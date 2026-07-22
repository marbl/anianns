"""Standalone NTRPrism k-mer spacing analysis and report generation."""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numba import njit, types
from numba.typed import Dict as NumbaDict

from anianns.kmer_pipeline import forward_kmer_hashes


@dataclass(frozen=True)
class SpacingPeak:
    """A group of nearby k-mer spacings represented as one histogram peak."""

    spacing: int
    minimum: int
    maximum: int
    count: int


@njit(cache=True)
def _valid_hash_distances(hashes, valid_windows):
    """Return positional gaps between consecutive valid occurrences of a hash."""
    last_indices = NumbaDict.empty(
        key_type=types.int32,
        value_type=types.int64,
    )
    distances = np.empty(len(hashes), dtype=np.int64)
    distance_count = 0
    for index in range(len(hashes)):
        if not valid_windows[index]:
            continue
        value = hashes[index]
        if value in last_indices:
            distances[distance_count] = index - last_indices[value]
            distance_count += 1
        last_indices[value] = index
    return distances[:distance_count].copy()


def _valid_kmer_windows(sequence: str, kmer: int) -> np.ndarray:
    """Mark k-mer starts whose complete window contains only A/C/G/T."""
    encoded = np.frombuffer(sequence.upper().encode("ascii"), dtype=np.uint8)
    if len(encoded) < kmer:
        return np.empty(0, dtype=np.bool_)
    valid_bases = (
        (encoded == ord("A"))
        | (encoded == ord("C"))
        | (encoded == ord("G"))
        | (encoded == ord("T"))
    )
    invalid_prefix = np.empty(len(encoded) + 1, dtype=np.int64)
    invalid_prefix[0] = 0
    np.cumsum(~valid_bases, out=invalid_prefix[1:])
    return (invalid_prefix[kmer:] - invalid_prefix[:-kmer]) == 0


def merge_spacing_counts(
    spacing_counts: Mapping[int, int], tolerance: int = 1
) -> list[SpacingPeak]:
    """Merge values around local count maxima, then rank peaks by count.

    The strongest unassigned spacing becomes a peak center and absorbs every
    still-unassigned value within ``tolerance`` bases. This avoids transitive
    chains (for example, a dense run of 1, 2, 3, ... must not become one giant
    peak) while still stacking neighboring values such as 170 and 171.
    """
    if tolerance < 0:
        raise ValueError("tolerance must be non-negative")
    observed = sorted(
        (int(spacing), int(count))
        for spacing, count in spacing_counts.items()
        if int(spacing) > 0 and int(count) > 0
    )
    if not observed:
        return []

    spacings = [spacing for spacing, _ in observed]
    counts = [count for _, count in observed]
    assigned = [False] * len(observed)
    centers = sorted(
        range(len(observed)),
        key=lambda index: (-counts[index], spacings[index]),
    )
    peaks = []
    for center_index in centers:
        if assigned[center_index]:
            continue
        center = spacings[center_index]
        left = bisect_left(spacings, center - tolerance)
        right = bisect_right(spacings, center + tolerance)
        member_indices = [index for index in range(left, right) if not assigned[index]]
        total = sum(counts[index] for index in member_indices)
        weighted_total = sum(
            spacings[index] * counts[index] for index in member_indices
        )
        representative = int(round(weighted_total / total))
        for index in member_indices:
            assigned[index] = True
        peaks.append(
            SpacingPeak(
                spacing=representative,
                minimum=min(spacings[index] for index in member_indices),
                maximum=max(spacings[index] for index in member_indices),
                count=total,
            )
        )
    return sorted(peaks, key=lambda peak: (-peak.count, peak.spacing))


def analyze_kmer_spacings(
    sequence: str, kmer: int = 21, merge_distance: int = 1
) -> tuple[list[SpacingPeak], int]:
    """Calculate ranked, merged spacing peaks for repeated forward k-mers."""
    if kmer <= 0:
        raise ValueError("k-mer length must be positive")
    if merge_distance < 0:
        raise ValueError("merge distance must be non-negative")
    if len(sequence) < kmer:
        raise ValueError(
            f"region length ({len(sequence)}) is shorter than k-mer length ({kmer})"
        )

    hashes = forward_kmer_hashes(sequence, kmer)
    valid_windows = _valid_kmer_windows(sequence, kmer)
    distances = _valid_hash_distances(hashes, valid_windows)
    if not len(distances):
        return [], 0
    values, counts = np.unique(distances, return_counts=True)
    raw_counts = {
        int(value): int(count)
        for value, count in zip(values.tolist(), counts.tolist())
    }
    return merge_spacing_counts(raw_counts, merge_distance), int(len(distances))


def format_spacing_table(
    peaks: list[SpacingPeak], interval_length: int, top_n: int = 10
) -> str:
    """Format the strongest spacing peaks and their share of the interval."""
    lines = [
        "Top k-mer spacing distances",
        "Rank  Distance (bp)  Merged range (bp)  Count      Interval %",
    ]
    for rank, peak in enumerate(peaks[:top_n], 1):
        merged_range = (
            str(peak.minimum)
            if peak.minimum == peak.maximum
            else f"{peak.minimum}-{peak.maximum}"
        )
        percentage = 100.0 * peak.count / interval_length if interval_length else 0.0
        lines.append(
            f"{rank:>4}  {peak.spacing:>13,}  {merged_range:>17}  "
            f"{peak.count:>9,}  {percentage:>10.3f}%"
        )
    if not peaks:
        lines.append("No repeated valid k-mers were found in the requested region.")
    return "\n".join(lines)


def format_ascii_histogram(
    peaks: list[SpacingPeak],
    interval_length: int,
    top_n: int = 10,
    width: int = 48,
) -> str:
    """Format a horizontal ASCII histogram for the strongest spacing peaks."""
    if width <= 0:
        raise ValueError("histogram width must be positive")
    shown = peaks[:top_n]
    lines = [
        "ASCII histogram (bars scaled to strongest peak; labels are % of interval)"
    ]
    if not shown:
        lines.append("(no repeated valid k-mers)")
        return "\n".join(lines)

    strongest_count = max(peak.count for peak in shown)
    label_width = max(len(f"{peak.spacing:,} bp") for peak in shown)
    for peak in shown:
        bar_length = max(1, round(width * peak.count / strongest_count))
        percentage = 100.0 * peak.count / interval_length if interval_length else 0.0
        lines.append(
            f"{f'{peak.spacing:,} bp':>{label_width}} | "
            f"{'#' * bar_length:<{width}} {percentage:>7.3f}%"
        )
    return "\n".join(lines)


def write_spacing_report(
    output_path: str | Path,
    *,
    fasta_path: str | Path,
    sequence_id: str,
    start: int,
    end: int,
    kmer: int,
    merge_distance: int,
    peaks: list[SpacingPeak],
    total_distances: int,
    top_n: int = 10,
) -> Path:
    """Write the ranked peaks and terminal-style histogram to a text report."""
    output_path = Path(output_path)
    lines = [
        "AniAnn's NTRPrism k-mer spacing report",
        f"FASTA\t{fasta_path}",
        f"Sequence\t{sequence_id}",
        f"Region_0_based_half_open\t{start}\t{end}",
        f"Kmer_length\t{kmer}",
        f"Nearby_merge_distance_bp\t{merge_distance}",
        f"Total_repeated_kmer_spacings\t{total_distances}",
        "",
        "Rank\tPeak_spacing_bp\tMerged_range_bp\tCount\t"
        "Interval_percentage\tRepeated_spacing_fraction",
    ]
    for rank, peak in enumerate(peaks[:top_n], 1):
        merged_range = (
            str(peak.minimum)
            if peak.minimum == peak.maximum
            else f"{peak.minimum}-{peak.maximum}"
        )
        interval_length = end - start
        interval_percentage = (
            100.0 * peak.count / interval_length if interval_length else 0.0
        )
        repeated_fraction = peak.count / total_distances if total_distances else 0.0
        lines.append(
            f"{rank}\t{peak.spacing}\t{merged_range}\t{peak.count}\t"
            f"{interval_percentage:.6f}\t{repeated_fraction:.6f}"
        )
    if not peaks:
        lines.append("No repeated valid k-mers were found in the requested region.")
    lines.extend(
        [
            "",
            format_ascii_histogram(peaks, end - start, top_n=top_n),
        ]
    )
    output_path.write_text("\n".join(lines) + "\n")
    return output_path


def write_spacing_histogram(
    output_path: str | Path,
    *,
    peaks: list[SpacingPeak],
    sequence_id: str,
    start: int,
    end: int,
    kmer: int,
    max_peaks: int = 50,
) -> Path:
    """Save a histogram-style bar plot of the strongest merged peaks."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter

    output_path = Path(output_path)
    shown = sorted(peaks[:max_peaks], key=lambda peak: peak.spacing)
    figure, axis = plt.subplots(figsize=(12, 6), dpi=160)
    if shown:
        positions = np.arange(len(shown))
        labels = [str(peak.spacing) for peak in shown]
        interval_length = end - start
        proportions = [
            peak.count / interval_length if interval_length else 0.0 for peak in shown
        ]
        strongest = {peak.spacing for peak in peaks[:10]}
        colors = [
            "#F4A261" if peak.spacing in strongest else "#1F7A8C"
            for peak in shown
        ]
        axis.bar(
            positions,
            proportions,
            width=0.85,
            color=colors,
            edgecolor="#16324F",
        )
        axis.set_xticks(positions, labels, rotation=90, fontsize=7)
    else:
        axis.text(
            0.5,
            0.5,
            "No repeated valid k-mers found",
            transform=axis.transAxes,
            ha="center",
            va="center",
            fontsize=14,
        )
    figure.suptitle(
        f"NTRPrism k-mer spacing peaks: {sequence_id}:{start}-{end} (k={kmer})"
    )
    axis.set_xlabel("Merged spacing peak (bp); orange bars are the report's top 10")
    axis.set_ylabel("Peak count as proportion of interval")
    axis.yaxis.set_major_formatter(PercentFormatter(xmax=1.0))
    axis.grid(axis="y", alpha=0.25)
    axis.margins(y=0.12)
    figure.tight_layout(rect=(0, 0, 1, 0.94))
    figure.savefig(output_path)
    plt.close(figure)
    return output_path
