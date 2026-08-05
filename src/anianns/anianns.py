import argparse
from contextlib import redirect_stderr, redirect_stdout
import csv
from datetime import datetime
from tokenize import group
from anianns.const import ASCII_ART, DESCRIPTION, VERSION
from itertools import islice
import polars as pl
import sys
import json
import os
import pysam
import math
import time
import numpy as np
from alive_progress import alive_bar
from numba import config as numba_config
from numba import set_num_threads

from anianns.ani_matrix import (
    intersection_matrix,
    intersection_matrix_rectangular,
    intersection_matrix_thresholded,
    intersection_matrix_with_threshold,
)

from anianns.build_kmer_db import (
    save_kmer_sets_shared_k,
)

from anianns.general_utils import (
    add_prefix_to_tuples,
    convert_dataframe_format,
    check_bed_vs_indexed_fasta,
    define_bounds,
    extract_region,
    extract_regions_by_name,
    extract_histograms_by_name,
    get_input_headers,
    plot_matrix,
    read_bed_files,
    validate_json,
    write_summary_file,
)

from anianns.kmer_utils import (
    build_kmer_sets,
    build_kmer_sets_multi,
)

from anianns.kmer_pipeline import (
    SequenceBandPlan,
    iter_hashed_fasta_bands,
    load_cached_sequence_hashes,
)

from anianns.ntrprism import (
    analyze_kmer_spacings,
    format_ascii_histogram,
    format_spacing_table,
    supports_periodic_lag,
    write_spacing_histogram,
    write_spacing_report,
)

from anianns.multi_window import (
    WindowCandidate,
    candidate_passes_support,
    derive_window_sizes,
    select_multi_window_candidates,
    write_window_selection,
)

from anianns.parse_matrix import (
    CandidateNeighborhoodAccumulator,
    DistalSatelliteLink,
    append_coordinates,
    deduplicate_distal_links,
    detect_adjacent_band_bridge,
    detect_candidate_to_all_links_from_full_matrix,
    detect_distal_links,
    filter_candidate_distal_links,
    get_diagonal_span,
    get_diagonal_span_from_sets,
    merge_shared_boundaries,
    PeriodicLagCandidate,
    detect_periodic_lag_candidates_from_matrix,
    detect_periodic_lag_candidates_from_sets,
    snap_distal_links_to_candidates,
    sobel_edge_response,
    sobel_with_diagonal_probes,
    sobel_spans,
)

from anianns.refine_boundaries import detect_precise_boundaries, report_borders

from anianns.union_find import (
    SatelliteDSU,
    sobel,
    sobel_with_diagonal_probes2,
    find_offdiag_rectangles,
)


MAX_THREADS = int(numba_config.NUMBA_NUM_THREADS)
ANNOTATE_LOG_FILENAME_PATTERN = "anianns_annotation_log_YYYY-MM-DD_HH-MM-SS.txt"


def annotate_log_filename(run_time=None):
    """Return a filesystem-safe log name containing the local run time."""
    run_time = datetime.now() if run_time is None else run_time
    return f"anianns_annotation_log_{run_time:%Y-%m-%d_%H-%M-%S}.txt"


class TeeTextStream:
    """Write text to both the terminal stream and a persistent log."""

    def __init__(self, terminal_stream, log_stream):
        self.terminal_stream = terminal_stream
        self.log_stream = log_stream

    def write(self, text):
        self.terminal_stream.write(text)
        self.log_stream.write(text)
        return len(text)

    def flush(self):
        self.terminal_stream.flush()
        self.log_stream.flush()

    def isatty(self):
        # Keep progress output readable in the plain-text log.
        return False

    def fileno(self):
        return self.terminal_stream.fileno()

    @property
    def encoding(self):
        return getattr(self.terminal_stream, "encoding", "utf-8")


def thread_count_type(value):
    """Parse a positive thread count supported by this Numba runtime."""
    try:
        threads = int(value)
    except (TypeError, ValueError) as error:
        raise argparse.ArgumentTypeError("thread count must be an integer") from error
    if threads < 1:
        raise argparse.ArgumentTypeError("thread count must be at least 1")
    if threads > MAX_THREADS:
        raise argparse.ArgumentTypeError(
            f"thread count cannot exceed the available maximum ({MAX_THREADS})"
        )
    return threads


def annotation_thread_allocation(threads, *, cache_hit, band_count):
    """Split the total budget between matrix work and streamed hashing."""
    use_hash_process = not cache_hit and threads > 1 and band_count > 1
    matrix_threads = threads - 1 if use_hash_process else threads
    return matrix_threads, use_hash_process


os.environ["KMP_WARNINGS"] = "FALSE"


def mask_type(value):
    """Allow either integer or string values for --mask."""
    try:
        return int(value)
    except ValueError:
        return value


def validate_ntrprism_range(start, end, seq_len):
    """
    Validate the --range argument for the ntrprism subcommand.

    Returns an error message string on failure, or None on success.
    """
    if start >= end:
        return f"[ERROR] --range start ({start}) must be less than end ({end})."
    if start < 0 or end > seq_len:
        return f"[ERROR] --range ({start}, {end}) is out of bounds for sequence of length {seq_len}."
    return None


def run_ntrprism_command(args):
    """Validate and analyze one region, printing results and optionally saving."""
    fasta_path = args.fasta
    if not os.path.isfile(fasta_path):
        print(f"[ERROR] FASTA file does not exist: {fasta_path}", file=sys.stderr)
        raise SystemExit(1)

    fasta = None
    try:
        try:
            fasta = pysam.FastaFile(fasta_path)
        except (OSError, ValueError) as error:
            print(
                f"[ERROR] Unable to open FASTA file {fasta_path}: {error}",
                file=sys.stderr,
            )
            raise SystemExit(1)

        if args.seq_id not in fasta.references:
            print(
                f"[ERROR] Sequence ID '{args.seq_id}' was not found in {fasta_path}.",
                file=sys.stderr,
            )
            raise SystemExit(1)

        sequence_length = fasta.get_reference_length(args.seq_id)
        start, end = args.range
        range_error = validate_ntrprism_range(start, end, sequence_length)
        if range_error:
            print(range_error, file=sys.stderr)
            raise SystemExit(1)

        try:
            sequence = fasta.fetch(args.seq_id, start, end)
        except (KeyError, OSError, ValueError) as error:
            print(
                f"[ERROR] Unable to fetch region {args.seq_id}:{start}-{end} "
                f"from {fasta_path}: {error}",
                file=sys.stderr,
            )
            raise SystemExit(1)
    finally:
        if fasta is not None:
            fasta.close()

    try:
        peaks, total_distances = analyze_kmer_spacings(
            sequence,
            kmer=args.kmer,
            merge_distance=args.merge_distance,
        )
    except ValueError as error:
        print(f"[ERROR] Cannot analyze requested region: {error}.", file=sys.stderr)
        raise SystemExit(1)

    interval_length = end - start
    if not args.quiet:
        print(
            f"NTRPrism: {args.seq_id}:{start}-{end} "
            f"({interval_length:,} bp, k={args.kmer})\n"
        )
        print(format_spacing_table(peaks, interval_length, top_n=10))
        print()
        print(format_ascii_histogram(peaks, interval_length, top_n=10))

    if not args.save:
        return None

    output_directory = args.directory or os.getcwd()
    try:
        os.makedirs(output_directory, exist_ok=True)
    except OSError as error:
        print(
            f"[ERROR] Unable to create output directory {output_directory}: {error}",
            file=sys.stderr,
        )
        raise SystemExit(1)

    safe_sequence_id = "".join(
        character if character.isalnum() or character in "._-" else "_"
        for character in args.seq_id
    )
    output_stem = f"{safe_sequence_id}_{start}_{end}_ntrprism"
    report_path = os.path.join(output_directory, f"{output_stem}.txt")
    histogram_path = os.path.join(output_directory, f"{output_stem}.png")
    try:
        write_spacing_report(
            report_path,
            fasta_path=fasta_path,
            sequence_id=args.seq_id,
            start=start,
            end=end,
            kmer=args.kmer,
            merge_distance=args.merge_distance,
            peaks=peaks,
            total_distances=total_distances,
            top_n=10,
        )
        write_spacing_histogram(
            histogram_path,
            peaks=peaks,
            sequence_id=args.seq_id,
            start=start,
            end=end,
            kmer=args.kmer,
        )
    except OSError as error:
        print(f"[ERROR] Unable to write NTRPrism output: {error}", file=sys.stderr)
        raise SystemExit(1)

    if not args.quiet:
        print(f"Saved NTRPrism top-10 spacing report to {report_path}")
        print(f"Saved NTRPrism spacing histogram to {histogram_path}")
    return report_path, histogram_path


def format_sequence_size(length):
    """Format a sequence length compactly with decimal genomic units."""
    length = int(length)
    if length < 0:
        raise ValueError("sequence length cannot be negative")
    if length >= 1_000_000_000:
        value, unit = length / 1_000_000_000, "gb"
    elif length >= 1_000_000:
        value, unit = length / 1_000_000, "mb"
    else:
        value, unit = length / 1000, "kb"

    if value >= 100:
        precision = 0
    elif value >= 10:
        precision = 1
    elif value >= 1:
        precision = 2
    else:
        precision = 3
    formatted = f"{value:.{precision}f}"
    if "." in formatted:
        formatted = formatted.rstrip("0").rstrip(".")
    return f"{formatted}{unit}"


def announce_matrix_creation(seq_id, sequence_length, cache_hit_dir=None):
    """Report a verified hash-cache hit immediately before matrix creation."""
    if cache_hit_dir is not None:
        print(f"Found hashes in {os.path.abspath(cache_hit_dir)}")
    print(
        f"Creating an ANI matrix for {seq_id} "
        f"({format_sequence_size(sequence_length)}):\n"
    )


def iter_with_stage_timing(iterator, timings, stage):
    """Accumulate only the time spent waiting for an iterator's next item."""
    while True:
        started = time.perf_counter()
        try:
            item = next(iterator)
        except StopIteration:
            timings[stage] += time.perf_counter() - started
            return
        timings[stage] += time.perf_counter() - started
        yield item


def scan_additional_window_candidates(
    hashed_band,
    band_plan,
    windows,
    *,
    kmer,
    identity,
    sketch,
    prebuilt_sketches=None,
):
    """Scan extra resolutions from shared band hashes without dense matrices."""
    if not windows:
        return [], 0.0, 0.0

    if prebuilt_sketches is None:
        sketches, sketch_runtime = build_band_window_sketches(
            hashed_band, band_plan, windows, sketch=sketch
        )
    else:
        sketches = prebuilt_sketches
        sketch_runtime = 0.0

    scan_runtime = 0.0
    candidates = []
    for window in windows:
        overlapping, non_overlapping = sketches[window]
        scan_started = time.perf_counter()
        spans = get_diagonal_span_from_sets(
            overlapping,
            non_overlapping,
            window,
            kmer,
            identity,
            zero_tol=2,
        )
        local_candidates = merge_shared_boundaries(
            intervals=spans,
            prefix=0,
            window=window,
            verbose=False,
        )
        scan_runtime += time.perf_counter() - scan_started
        candidates.extend(
            WindowCandidate(
                start=start + hashed_band.start,
                end=end + hashed_band.start,
                count=count,
                window=window,
            )
            for start, end, count in local_candidates
        )
    return candidates, sketch_runtime, scan_runtime


def scan_primary_window(
    overlapping,
    non_overlapping,
    window,
    kmer,
    identity,
    *,
    threshold_matrix=None,
):
    """Find primary diagonal and periodic calls, reusing a dense matrix."""
    if threshold_matrix is None:
        spans = get_diagonal_span_from_sets(
            overlapping,
            non_overlapping,
            window,
            kmer,
            identity,
            zero_tol=2,
        )
        periodic = detect_periodic_lag_candidates_from_sets(
            overlapping,
            non_overlapping,
            window,
            kmer,
            identity,
        )
    else:
        spans = get_diagonal_span(threshold_matrix, window, zero_tol=2)
        periodic = detect_periodic_lag_candidates_from_matrix(
            threshold_matrix,
            window,
        )
    return spans, periodic


def combine_band_candidates(primary_candidates, additional_candidates, prefix):
    """Combine primary and supported multi-window calls in band coordinates."""
    combined = [tuple(candidate) for candidate in primary_candidates]
    seen = {(int(start), int(end)) for start, end, *_rest in combined}
    for candidate in additional_candidates:
        if not candidate_passes_support(candidate):
            continue
        local = (
            int(candidate.start) - int(prefix),
            int(candidate.end) - int(prefix),
            int(candidate.count),
        )
        key = local[:2]
        if key not in seen:
            combined.append(local)
            seen.add(key)
    return combined


def merge_periodic_lag_candidates(candidates, window):
    """Merge compatible periodic rescues split at matrix-band boundaries."""
    merged = []
    for candidate in sorted(candidates, key=lambda item: (item.start, item.end)):
        if not merged:
            merged.append(candidate)
            continue
        previous = merged[-1]
        period_difference = abs(previous.period_bp - candidate.period_bp)
        period_tolerance = max(window, int(0.10 * min(previous.period_bp, candidate.period_bp)))
        if candidate.start <= previous.end + (2 * window) and period_difference <= period_tolerance:
            total_length = (previous.end - previous.start) + (
                candidate.end - candidate.start
            )
            coverage = (
                (previous.coverage * (previous.end - previous.start))
                + (candidate.coverage * (candidate.end - candidate.start))
            ) / max(1, total_length)
            stronger = max(
                (previous, candidate),
                key=lambda item: item.coverage * item.contrast,
            )
            merged[-1] = PeriodicLagCandidate(
                start=previous.start,
                end=max(previous.end, candidate.end),
                count=previous.count + candidate.count,
                period_bp=stronger.period_bp,
                harmonic_lags=stronger.harmonic_lags,
                coverage=coverage,
                contrast=min(previous.contrast, candidate.contrast),
            )
        else:
            merged.append(candidate)
    return merged


def _covered_fraction(candidate, existing_candidates):
    overlaps = []
    for existing in existing_candidates:
        start = max(candidate.start, existing.start)
        end = min(candidate.end, existing.end)
        if end > start:
            overlaps.append((start, end))
    covered = 0
    if overlaps:
        current_start, current_end = sorted(overlaps)[0]
        for start, end in sorted(overlaps)[1:]:
            if start > current_end:
                covered += current_end - current_start
                current_start, current_end = start, end
            else:
                current_end = max(current_end, end)
        covered += current_end - current_start
    return covered / max(1, candidate.end - candidate.start)


def validate_periodic_lag_candidates(
    fasta_handle,
    sequence_id,
    proposals,
    existing_candidates,
    *,
    window,
    kmer,
    verbose=False,
):
    """Validate uncovered matrix-stripe proposals with sequence spacing peaks."""
    validated = []
    for proposal in merge_periodic_lag_candidates(proposals, window):
        if _covered_fraction(proposal, existing_candidates) >= 0.80:
            continue
        sequence = fasta_handle.fetch(sequence_id, proposal.start, proposal.end)
        if not supports_periodic_lag(
            sequence,
            proposal.period_bp,
            kmer=kmer,
        ):
            continue
        validated.append(
            WindowCandidate(
                proposal.start,
                proposal.end,
                proposal.count,
                window,
                source="periodic_diagonal",
            )
        )
        if verbose:
            print(
                "Periodic diagonal rescue: "
                f"{sequence_id}:{proposal.start:,}-{proposal.end:,} "
                f"period≈{proposal.period_bp:,} bp "
                f"coverage={proposal.coverage:.1%}"
            )
    return validated


def candidate_genomic_ranges(candidates, prefix):
    """Convert band-local satellite candidates to genomic intervals."""
    return [
        (int(prefix) + int(start), int(prefix) + int(end))
        for start, end, *_rest in candidates
        if int(end) > int(start)
    ]


def build_band_window_sketches(hashed_band, band_plan, windows, *, sketch):
    """Build all requested band resolutions with one modulo-mask pass."""
    window_configs = {}
    for window in windows:
        window_plan = band_plan.window_plans[window]
        hash_count = len(band_plan.hashes_for_window(hashed_band, window))
        window_configs[window] = (
            hash_count,
            window_plan.max_len,
            window_plan.interval,
        )

    started = time.perf_counter()
    sketches = build_kmer_sets_multi(
        hashed_band.hashes,
        window_configs,
        sketch=sketch,
    )
    return sketches, time.perf_counter() - started


def save_matrix_heatmap(
    matrix,
    output_directory,
    seq_id,
    matrix_index,
    genomic_start,
    window,
    identity,
    distal_links=(),
    satellite_ranges=(),
    identity_matrix=None,
):
    """Save annotated Sobel and high-resolution identity views of one band."""
    plot_directory = os.path.join(output_directory, "matrix_plots")
    os.makedirs(plot_directory, exist_ok=True)
    safe_seq_id = str(seq_id).replace(os.sep, "_")
    if os.altsep:
        safe_seq_id = safe_seq_id.replace(os.altsep, "_")
    plot_path = os.path.join(
        plot_directory,
        f"{safe_seq_id}_matrix_{matrix_index:04d}.png",
    )
    identity_path = os.path.join(
        plot_directory,
        f"{safe_seq_id}_matrix_{matrix_index:04d}_identity.png",
    )
    genomic_end = genomic_start + (matrix.shape[0] * window)
    display_matrix, display_scale = downsample_binary_matrix(matrix)
    edge_overlay = sobel_edge_response(display_matrix)
    diagonal_ranges = [
        (int(start) - genomic_start, int(end) - genomic_start)
        for start, end in satellite_ranges
        if int(end) > int(start)
    ]
    highlight_ranges = []
    for link in distal_links:
        # The similarity matrix is symmetric, so show both orientations.
        highlight_ranges.extend(
            (
                (
                    link.start2 - genomic_start,
                    link.end2 - genomic_start,
                    link.start1 - genomic_start,
                    link.end1 - genomic_start,
                ),
                (
                    link.start1 - genomic_start,
                    link.end1 - genomic_start,
                    link.start2 - genomic_start,
                    link.end2 - genomic_start,
                ),
            )
        )
    plot_matrix(
        display_matrix,
        title=(
            f"{seq_id} matrix {matrix_index}: "
            f"{genomic_start:,}-{genomic_end:,} bp"
        ),
        show_colorbar=False,
        dpi=100,
        figsize=(5.5, 5.5),
        save_path=plot_path,
        offset=window * display_scale,
        coordinate_origin=genomic_start,
        coordinate_end=genomic_end,
        coordinate_units="Mbp",
        vmin=0,
        vmax=1,
        aspect="equal",
        highlight_ranges=highlight_ranges,
        diagonal_ranges=diagonal_ranges,
        edge_overlay=edge_overlay,
        edge_only=True,
        legend_outside=True,
    )
    if identity_matrix is None:
        identity_matrix = np.where(matrix, 100.0, 0.0)
    identity_display, identity_scale = downsample_numeric_matrix(
        identity_matrix, max_pixels=1024
    )
    plot_matrix(
        identity_display,
        title=(
            f"{seq_id} matrix {matrix_index}: "
            f"{genomic_start:,}-{genomic_end:,} bp"
        ),
        cmap="spectral_11_r",
        show_colorbar=True,
        colorbar_label="ANI identity (%)",
        dpi=160,
        figsize=(8.5, 7.0),
        save_path=identity_path,
        offset=window * identity_scale,
        coordinate_origin=genomic_start,
        coordinate_end=genomic_end,
        coordinate_units="Mbp",
        colorbar_pad=0.12,
        reserve_colorbar_space=True,
        vmin=identity,
        vmax=100,
        white_below=identity,
        aspect="equal",
    )
    return plot_path


def save_adjacent_matrix_heatmap(
    previous_matrix,
    current_matrix,
    cross_matrix,
    output_directory,
    seq_id,
    previous_index,
    current_index,
    genomic_start,
    window,
    identity,
    distal_links=(),
    seam_candidates=(),
    satellite_ranges=(),
    previous_identity_matrix=None,
    current_identity_matrix=None,
    cross_identity_matrix=None,
):
    """Save annotated Sobel and identity views of an adjacent band pair."""
    pair_directory = os.path.join(output_directory, "matrix_pairs")
    os.makedirs(pair_directory, exist_ok=True)
    safe_seq_id = str(seq_id).replace(os.sep, "_")
    if os.altsep:
        safe_seq_id = safe_seq_id.replace(os.altsep, "_")
    plot_path = os.path.join(
        pair_directory,
        f"{safe_seq_id}_matrices_{previous_index:04d}_{current_index:04d}.png",
    )
    identity_path = os.path.join(
        pair_directory,
        f"{safe_seq_id}_matrices_{previous_index:04d}_{current_index:04d}_identity.png",
    )
    combined = np.block(
        [
            [previous_matrix, cross_matrix],
            [cross_matrix.T, current_matrix],
        ]
    )
    display_matrix, display_scale = downsample_binary_matrix(combined)
    edge_overlay = sobel_edge_response(display_matrix)
    genomic_end = genomic_start + (combined.shape[0] * window)
    diagonal_ranges = [
        (int(start) - genomic_start, int(end) - genomic_start)
        for start, end in satellite_ranges
        if int(end) > int(start)
    ]
    highlight_ranges = []
    for link in distal_links:
        highlight_ranges.extend(
            (
                (
                    link.start2 - genomic_start,
                    link.end2 - genomic_start,
                    link.start1 - genomic_start,
                    link.end1 - genomic_start,
                ),
                (
                    link.start1 - genomic_start,
                    link.end1 - genomic_start,
                    link.start2 - genomic_start,
                    link.end2 - genomic_start,
                ),
            )
        )
    for start, end, _count in seam_candidates:
        diagonal_ranges.append((start - genomic_start, end - genomic_start))
    plot_matrix(
        display_matrix,
        title=(
            f"{seq_id}\nmatrices {previous_index}-{current_index}: "
            f"{genomic_start:,}-{genomic_end:,} bp"
        ),
        show_colorbar=False,
        dpi=100,
        figsize=(5.5, 5.5),
        save_path=plot_path,
        offset=window * display_scale,
        coordinate_origin=genomic_start,
        coordinate_end=genomic_end,
        coordinate_units="Mbp",
        vmin=0,
        vmax=1,
        aspect="equal",
        highlight_ranges=highlight_ranges,
        diagonal_ranges=diagonal_ranges,
        edge_overlay=edge_overlay,
        edge_only=True,
        legend_outside=True,
    )
    if previous_identity_matrix is None:
        previous_identity_matrix = np.where(previous_matrix, 100.0, 0.0)
    if current_identity_matrix is None:
        current_identity_matrix = np.where(current_matrix, 100.0, 0.0)
    if cross_identity_matrix is None:
        cross_identity_matrix = np.where(cross_matrix, identity, 0.0)
    identity_combined = np.block(
        [
            [previous_identity_matrix, cross_identity_matrix],
            [cross_identity_matrix.T, current_identity_matrix],
        ]
    )
    identity_display, identity_scale = downsample_numeric_matrix(
        identity_combined, max_pixels=1024
    )
    plot_matrix(
        identity_display,
        title=(
            f"{seq_id}\nmatrices {previous_index}-{current_index}: "
            f"{genomic_start:,}-{genomic_end:,} bp"
        ),
        cmap="spectral_11_r",
        show_colorbar=True,
        colorbar_label="ANI identity (%)",
        dpi=160,
        figsize=(8.5, 7.0),
        save_path=identity_path,
        offset=window * identity_scale,
        coordinate_origin=genomic_start,
        coordinate_end=genomic_end,
        coordinate_units="Mbp",
        colorbar_pad=0.12,
        reserve_colorbar_space=True,
        vmin=identity,
        vmax=100,
        white_below=identity,
        aspect="equal",
    )
    return plot_path


def incorporate_seam_candidates(candidates, seam_candidates, window):
    """Merge boundary-spanning evidence with candidate fragments in place."""
    for seam_start, seam_end, seam_count in seam_candidates:
        matching = []
        for index, (start, end, _count) in enumerate(candidates):
            gap = max(start - seam_end, seam_start - end, 0)
            if gap <= 2 * window and start < seam_end and end > seam_start:
                matching.append(index)
        if matching:
            matched = [candidates[index] for index in matching]
            merged_start = min(seam_start, *(item[0] for item in matched))
            merged_end = max(seam_end, *(item[1] for item in matched))
            merged_count = max(
                seam_count,
                sum(item[2] for item in matched),
            )
            for index in reversed(matching):
                del candidates[index]
            candidates.append((merged_start, merged_end, merged_count))
        else:
            candidates.append((seam_start, seam_end, seam_count))
    candidates.sort(key=lambda candidate: (candidate[0], candidate[1]))


def detect_and_save_matrix_heatmap(
    matrix,
    output_directory,
    seq_id,
    matrix_index,
    genomic_start,
    window,
    identity,
    candidates,
):
    """Detect local distal blocks and overlay them on a saved heatmap."""
    links = detect_matrix_distal_links(
        matrix, genomic_start, window, candidates
    )
    plot_path = save_matrix_heatmap(
        matrix,
        output_directory,
        seq_id,
        matrix_index,
        genomic_start,
        window,
        identity,
        distal_links=links,
        satellite_ranges=candidate_genomic_ranges(candidates, genomic_start),
    )
    return plot_path, links


def detect_matrix_distal_links(matrix, genomic_start, window, candidates):
    """Detect and filter distal links in one band without rendering it."""
    genomic_starts = genomic_start + (
        np.arange(matrix.shape[0], dtype=np.int64) * window
    )
    candidate_ranges = tuple(
        (genomic_start + int(start), genomic_start + int(end))
        for start, end, _count in candidates
    )
    links = detect_distal_links(
        matrix,
        genomic_starts,
        candidate_ranges,
        window,
    )
    links = filter_candidate_distal_links(
        links,
        tuple((*interval, 0) for interval in candidate_ranges),
        proximity=2 * window,
    )
    return links


def downsample_binary_matrix(matrix, max_pixels=512):
    """Max-pool a boolean matrix to the plot's actual low-resolution raster."""
    binary = np.asarray(matrix, dtype=np.bool_)
    if binary.ndim != 2:
        raise ValueError("plot matrix must be two-dimensional")
    if max_pixels <= 0:
        raise ValueError("max_pixels must be positive")
    largest_axis = max(binary.shape, default=0)
    scale = max(1, math.ceil(largest_axis / max_pixels))
    if scale == 1 or binary.size == 0:
        return binary, scale

    padded_rows = math.ceil(binary.shape[0] / scale) * scale
    padded_columns = math.ceil(binary.shape[1] / scale) * scale
    padded = np.zeros((padded_rows, padded_columns), dtype=np.bool_)
    padded[: binary.shape[0], : binary.shape[1]] = binary
    pooled = padded.reshape(
        padded_rows // scale,
        scale,
        padded_columns // scale,
        scale,
    ).any(axis=(1, 3))
    return pooled, scale


def downsample_numeric_matrix(matrix, max_pixels=1024):
    """Max-pool a numeric identity matrix for a high-resolution plot raster."""
    values = np.asarray(matrix)
    if values.ndim != 2:
        raise ValueError("plot matrix must be two-dimensional")
    if max_pixels <= 0:
        raise ValueError("max_pixels must be positive")
    largest_axis = max(values.shape, default=0)
    scale = max(1, math.ceil(largest_axis / max_pixels))
    if scale == 1 or values.size == 0:
        return values, scale

    padded_rows = math.ceil(values.shape[0] / scale) * scale
    padded_columns = math.ceil(values.shape[1] / scale) * scale
    padded = np.zeros((padded_rows, padded_columns), dtype=values.dtype)
    padded[: values.shape[0], : values.shape[1]] = values
    pooled = padded.reshape(
        padded_rows // scale,
        scale,
        padded_columns // scale,
        scale,
    ).max(axis=(1, 3))
    return pooled, scale


def write_distal_links(links, seq_id, output_path, coordinate_offset=0):
    """Write detected distal relationships as BEDPE plus support metrics."""
    with open(output_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            (
                "#chrom1",
                "start1",
                "end1",
                "chrom2",
                "start2",
                "end2",
                "name",
                "score",
                "strand1",
                "strand2",
                "source",
                "density",
                "row_coverage",
                "column_coverage",
                "hit_count",
            )
        )
        for index, link in enumerate(links, start=1):
            writer.writerow(
                (
                    seq_id,
                    link.start1 + coordinate_offset,
                    link.end1 + coordinate_offset,
                    seq_id,
                    link.start2 + coordinate_offset,
                    link.end2 + coordinate_offset,
                    f"distal_link_{index:04d}",
                    min(1000, round(link.density * 1000)),
                    ".",
                    ".",
                    link.source,
                    f"{link.density:.6f}",
                    f"{link.row_coverage:.6f}",
                    f"{link.column_coverage:.6f}",
                    link.hit_count,
                )
            )


def promote_unmatched_distal_candidates(
    links,
    starts,
    ends,
    names,
    monomers,
    periodicities,
    hor_flags,
    *,
    fasta_file,
    seq_id,
    seq_len,
    window,
    k,
    sequence_hashes=None,
    verbose=False,
):
    """NTR-validate and add the missing endpoint of one-sided distal links."""

    def overlap(first, second):
        amount = max(0, min(first[1], second[1]) - max(first[0], second[0]))
        if amount == 0:
            return 0.0
        return amount / min(first[1] - first[0], second[1] - second[0])

    def match(interval):
        best_index = None
        best_overlap = 0.0
        for index, candidate in enumerate(zip(starts, ends)):
            candidate_overlap = overlap(interval, candidate)
            if candidate_overlap > best_overlap:
                best_index = index
                best_overlap = candidate_overlap
        return best_index if best_overlap >= 0.5 else None

    promoted_cache = {}
    validated_links = []
    ordered_links = sorted(
        links,
        key=lambda link: (
            link.source not in ("candidate", "candidate_to_all"),
            -(link.density + min(link.row_coverage, link.column_coverage)),
        ),
    )
    for link in ordered_links:
        intervals = [(link.start1, link.end1), (link.start2, link.end2)]
        matches = [match(interval) for interval in intervals]
        for axis, candidate_index in enumerate(matches):
            if candidate_index is not None:
                intervals[axis] = (
                    int(starts[candidate_index]),
                    int(ends[candidate_index]),
                )
        matched_count = sum(index is not None for index in matches)
        if matched_count == 0:
            continue
        if matched_count == 1:
            unknown_axis = 0 if matches[0] is None else 1
            unknown_interval = intervals[unknown_axis]
            if unknown_interval not in promoted_cache:
                previous = max(
                    (
                        (int(start), int(end))
                        for start, end in zip(starts, ends)
                        if int(end) <= unknown_interval[0]
                    ),
                    default=(0, 1),
                    key=lambda candidate: candidate[1],
                )
                try:
                    result = detect_precise_boundaries(
                        fasta_file=fasta_file,
                        seq_id=seq_id,
                        seq_len=seq_len,
                        window=window,
                        k=k,
                        coordinates=unknown_interval,
                        verbose=verbose,
                        classify=False,
                        previous_coordinates=previous,
                        sequence_hashes=sequence_hashes,
                    )
                except Exception as error:
                    if verbose:
                        print(
                            f"Unable to validate distal-only candidate "
                            f"{unknown_interval}: {error}"
                        )
                    result = None

                if result is None or result[3] is None or result[3][0] in (None, 0):
                    promoted_cache[unknown_interval] = None
                else:
                    promoted_start, promoted_end, promoted_name, prism = result
                    refined_interval = (int(promoted_start), int(promoted_end))
                    # Do not turn a fragment of the known endpoint into a new node.
                    if match(refined_interval) is not None:
                        promoted_cache[unknown_interval] = None
                    else:
                        starts.append(refined_interval[0])
                        ends.append(refined_interval[1])
                        names.append(promoted_name)
                        monomers.append(prism[0])
                        periodicities.append(prism[2])
                        hor_flags.append(bool(prism[1]))
                        promoted_cache[unknown_interval] = refined_interval
                        if verbose:
                            print(
                                f"Promoted distal-supported satellite at "
                                f"{refined_interval[0]}-{refined_interval[1]}"
                            )

            promoted_interval = promoted_cache.get(unknown_interval)
            if promoted_interval is None:
                continue
            intervals[unknown_axis] = promoted_interval

        first, second = intervals
        if first[0] > second[0]:
            first, second = second, first
        if match(first) is None or match(second) is None or match(first) == match(second):
            continue
        validated_links.append(
            DistalSatelliteLink(
                start1=first[0],
                end1=first[1],
                start2=second[0],
                end2=second[1],
                source=link.source,
                density=link.density,
                row_coverage=link.row_coverage,
                column_coverage=link.column_coverage,
                hit_count=link.hit_count,
            )
        )

    satellites = sorted(
        zip(starts, ends, names, monomers, periodicities, hor_flags),
        key=lambda satellite: (satellite[0], satellite[1]),
    )
    if satellites:
        (
            sorted_starts,
            sorted_ends,
            sorted_names,
            sorted_monomers,
            sorted_periodicities,
            sorted_hor_flags,
        ) = map(list, zip(*satellites))
        starts[:] = sorted_starts
        ends[:] = sorted_ends
        names[:] = sorted_names
        monomers[:] = sorted_monomers
        periodicities[:] = sorted_periodicities
        hor_flags[:] = sorted_hor_flags
    return deduplicate_distal_links(validated_links)


def resolve_overlapping_satellites(
    starts,
    ends,
    names,
    monomers,
    periodicities,
    hor_flags,
    *,
    verbose=False,
):
    """Make refined satellite intervals sorted and non-overlapping in place.

    Boundary refinement and distal promotion operate on each candidate
    independently. Two unrelated candidates can therefore acquire slightly
    conflicting boundary estimates even though neither is a duplicate of the
    other. Split a partial overlap halfway between those estimates, which makes
    the smallest symmetric adjustment. If one call is wholly contained in an
    earlier, longer call, discard the contained call because it cannot be
    represented as a separate non-overlapping interval.
    """
    values = (starts, ends, names, monomers, periodicities, hor_flags)
    lengths = {len(value) for value in values}
    if len(lengths) != 1:
        raise ValueError("satellite metadata lists must have equal lengths")

    satellites = sorted(
        zip(starts, ends, names, monomers, periodicities, hor_flags),
        key=lambda satellite: (
            int(satellite[0]),
            -(int(satellite[1]) - int(satellite[0])),
        ),
    )
    resolved = []
    for satellite in satellites:
        start, end, name, monomer, periodicity, is_hor = satellite
        current = [int(start), int(end), name, monomer, periodicity, is_hor]
        if current[1] <= current[0]:
            continue
        if not resolved or current[0] >= resolved[-1][1]:
            resolved.append(current)
            continue

        previous = resolved[-1]
        if current[1] <= previous[1]:
            if verbose:
                print(
                    "Discarding contained refined satellite "
                    f"{current[0]}-{current[1]} inside "
                    f"{previous[0]}-{previous[1]}."
                )
            continue

        old_previous_end = previous[1]
        old_current_start = current[0]
        boundary = (old_previous_end + old_current_start) // 2
        previous[1] = boundary
        current[0] = boundary
        if verbose:
            print(
                "Resolved overlapping satellite boundaries "
                f"{previous[0]}-{old_previous_end} and "
                f"{old_current_start}-{current[1]} at {boundary}."
            )
        resolved.append(current)

    columns = list(zip(*resolved)) if resolved else [[] for _ in range(6)]
    for destination, column in zip(values, columns):
        destination[:] = list(column)
    return tuple(values)


def get_parser():
    """
    Argument parsing for stand-alone runs.

    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=DESCRIPTION,
    )
    subparsers = parser.add_subparsers(
        dest="command", help="Choose mode: annotate, build_db, ntrprism"
    )
    annotate_parser = subparsers.add_parser(
        "annotate",
        help="Takes input fasta(s) and outputs an annotated bedfile of satellite arrays.",
    )
    build_db_parser = subparsers.add_parser(
        "build_db",
        help="Takes input fasta(s), bedfile(s) of known satellite coordinates, and a config file, and outputs a kme db (directory).",
    )
    ntrprism_parser = subparsers.add_parser(
        "ntrprism",
        help="Takes input fasta(s) plus regions of interest, outputs ntrprism k-mer spectra.",
    )
    annotate_parser.add_argument(
        "-f",
        "--fasta",
        default=argparse.SUPPRESS,
        help="Path to input fasta file(s).",
        required=True,
        nargs="+",
    )
    annotate_parser.add_argument(
        "-s",
        "--seq_id",
        nargs="+",
        default=None,
        help="Sequence ID to extract (multiple if using multifasta file). Will ignore if not found.",
    )
    annotate_parser.add_argument(
        "-d",
        "--directory",
        default=None,
        help="Name of output directory. Default: current working directory.",
    )
    annotate_parser.add_argument(
        "-o",
        "--output-format",
        choices=["bed", "gtf", "gff", "csv", "tsv", "json"],
        default="bed",
        help=(
            "Specify the output file format. "
            "Accepted values: bed, gtf, gff, csv, tsv, or json. "
            "Default: bed."
        ),
    )
    annotate_parser.add_argument(
        "-m",
        "--mask",
        nargs="*",
        type=mask_type,
        default=None,
        help=(
            "Name(s) or ID(s) of satellite arrays to mask. "
            "Use without parameters to apply default masking. "
            "Use 'ALL' for everything. Default: not applied."
        ),
    )
    annotate_parser.add_argument(
        "--soft",
        default=False,
        action="store_true",
        help="Apply soft masking (requires --mask).",
    )
    annotate_parser.add_argument(
        "-c",
        "--classify",
        default=None,
        help="Directory of satellite kmer db files, which Ani Ann's will use to classify into known satellite classes.",
    )
    annotate_parser.add_argument(
        "-t",
        "--threshold",
        default=None,
        help="Directory of satellite kmer db files, which Ani Ann's will use to classify into known satellite classes.",
    )
    annotate_parser.add_argument(
        "-k", "--kmer", type=int, default=21, help="k-mer length. Default: 21"
    )
    annotate_parser.add_argument(
        "--sketch",
        type=int,
        choices=(2, 4),
        default=4,
        help=(
            "Keep hashes divisible by this modulo when building window sketches. "
            "Use 2 for a denser, slower sketch or 4 for the default."
        ),
    )
    annotate_parser.add_argument(
        "-i", "--identity", type=int, default=86, help="Identity threshold. Default: 86"
    )
    annotate_parser.add_argument(
        "-w",
        "--window",
        type=int,
        default=2000,
        help=(
            "Central dotplot window size. AniAnn's also scans half and double "
            "this value using shared hashes. Default: 2000 (scans 1000, "
            "2000, and 4000)."
        ),
    )
    annotate_parser.add_argument(
        "--band",
        type=float,
        default=2.0,
        help="Max height in Mbp of band. Default: 2.0",
    )
    annotate_parser.add_argument(
        "--cache-dir",
        default=None,
        help=(
            "Opt-in directory for reusable canonical k-mer hash caches. "
            "A matching cache is loaded when present and created when absent. "
            "Default: caching disabled."
        ),
    )
    annotate_parser.add_argument(
        "-j",
        "--threads",
        type=thread_count_type,
        default=MAX_THREADS,
        help="Maximum total compute threads. Default: all available threads.",
    )
    annotate_parser.add_argument(
        "--identifier",
        help="Name of identifier. Used when no matches to a k-mer db are found, or if `--classify` is not provided. bed file to output to. Default: None",
    )
    annotate_parser.add_argument(
        "-p",
        "--plot",
        action="store_true",
        help=(
            "Save annotated Sobel and high-resolution Spectral identity views "
            "for each band under <output directory>/matrix_plots."
        ),
    )
    annotate_parser.add_argument(
        "--distal",
        action="store_true",
        help=(
            "Detect distal satellite links, validate unmatched axes, and union "
            "linked satellites in the DSU. With --plot, detected links are "
            "outlined in red."
        ),
    )
    annotate_group = annotate_parser.add_mutually_exclusive_group()
    annotate_group.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose logging output."
    )
    annotate_group.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress all logging output and text.",
    )
    annotate_parser.add_argument(
        "--log",
        action="store_true",
        help=(
            "Write verbose output to a timestamped file such as "
            f"{ANNOTATE_LOG_FILENAME_PATTERN} in the output directory. "
            "Implies --verbose."
        ),
    )

    build_db_parser.add_argument(
        "-f",
        "--fasta",
        default=argparse.SUPPRESS,
        required=True,
        help="Path to input fasta file(s).",
        nargs="+",
    )
    build_db_parser.add_argument(
        "-b",
        "--bed",
        default=argparse.SUPPRESS,
        required=True,
        help="Path to input bed file(s).",
        nargs="+",
    )
    build_db_parser.add_argument(
        "-c",
        "--config",
        help="Path to input config file(s).",
    )
    build_db_parser.add_argument(
        "-d",
        "--directory",
        default=None,
        help="Name of kmer db output directory. Default: In current working directory, will create a directory with the k-mer length.",
    )
    build_db_parser.add_argument(
        "-k", "--kmer", type=int, default=21, help="k-mer length. Default: 21"
    )
    build_db_parser.add_argument(
        "-j",
        "--threads",
        type=thread_count_type,
        default=MAX_THREADS,
        help="Maximum compute threads. Default: all available threads.",
    )
    build_db_group = build_db_parser.add_mutually_exclusive_group()
    build_db_group.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose logging output."
    )
    build_db_group.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress all logging output and text.",
    )

    ntrprism_parser.add_argument(
        "-f",
        "--fasta",
        help="Path to one indexed FASTA file.",
        required=True,
    )
    ntrprism_parser.add_argument(
        "-s",
        "--seq_id",
        required=True,
        help="Sequence ID containing the region of interest.",
    )
    ntrprism_parser.add_argument(
        "-r",
        "--range",
        required=True,
        type=int,
        nargs=2,
        metavar=("START", "END"),
        help="Region as 0-based, half-open coordinates: START END.",
    )
    ntrprism_parser.add_argument(
        "-d",
        "--directory",
        default=None,
        help=(
            "Output directory used with --save. Default: current directory."
        ),
    )
    ntrprism_parser.add_argument(
        "-k",
        "--kmer",
        type=int,
        default=21,
        help="K-mer length used to calculate repeated-k-mer spacings. Default: 21.",
    )
    ntrprism_parser.add_argument(
        "-j",
        "--threads",
        type=thread_count_type,
        default=MAX_THREADS,
        help="Maximum compute threads. Default: all available threads.",
    )
    ntrprism_parser.add_argument(
        "--merge-distance",
        type=int,
        default=1,
        help="Maximum gap in bp between neighboring spacing values merged into one peak.",
    )
    ntrprism_parser.add_argument(
        "--save",
        action="store_true",
        help="Save the text report and PNG histogram in addition to terminal output.",
    )
    ntrprism_parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress the terminal report (use with --save for files only).",
    )

    return parser


def main():
    args = get_parser().parse_args()
    if getattr(args, "quiet", False) and getattr(args, "log", False):
        return 2
    if getattr(args, "quiet", False):
        with open(os.devnull, "w") as null_stream:
            with redirect_stdout(null_stream), redirect_stderr(null_stream):
                return _run_command(args)

    if getattr(args, "log", False):
        args.verbose = True
        log_directory = args.directory or os.getcwd()
        os.makedirs(log_directory, exist_ok=True)
        log_path = os.path.join(log_directory, annotate_log_filename())
        with open(log_path, "w") as log_stream:
            with redirect_stdout(TeeTextStream(sys.stdout, log_stream)):
                with redirect_stderr(TeeTextStream(sys.stderr, log_stream)):
                    return _run_visible_command(args)

    return _run_visible_command(args)


def _run_visible_command(args):
    print(ASCII_ART)
    print(f" {VERSION}")
    print("─" * 65)
    return _run_command(args)


def _run_command(args):
    if hasattr(args, "threads"):
        set_num_threads(args.threads)

    # -------- BUILD DB LOGIC --------#
    if args.command == "build_db":
        if not args.directory:
            foldername = f"k_{args.kmer}"
            args.directory = os.path.join(os.getcwd(), foldername)
        if not args.quiet:
            print(f"Building a k-mer database....\n")
        if args.config:
            with open(args.config, "r") as f:
                if not validate_json(args.config):
                    sys.exit(1)
                else:
                    try:
                        satellite_metadata = json.load(f)
                        print(f"Successfully loaded sat metadata from {args.config}!\n")
                    except:
                        print(
                            f"Error loading satellite metadata from  {args.config}. Exiting.\n"
                        )
                        sys.exit(1)
        else:
            satellite_metadata = None
            if not args.quiet:
                print(
                    f"No metadata provided. Creating a compressed kmer db for each unique satellite.\n"
                )

        if not args.quiet:
            print(f"Creating kmer db using {args.fasta}\n")
        bedfiles = read_bed_files(args.bed)

        check_bed_vs_indexed_fasta(bedfiles, args.fasta)

        # Normalize all BED DataFrames to have consistent schema before concatenating
        normalized_bedfiles = []
        required_columns = ["chrom", "start", "end", "name"]

        for df in bedfiles:
            # Ensure all required columns exist
            for col in required_columns:
                if col not in df.columns:
                    if col == "name":
                        # Use chrom as default name if missing
                        df = df.with_columns(pl.col("chrom").alias("name"))
                    else:
                        raise ValueError(
                            f"Required column '{col}' missing from BED file"
                        )

            # Select only the columns we need for consistency
            normalized_df = df.select(required_columns)
            normalized_bedfiles.append(normalized_df)

        # Concatenate normalized DataFrames
        bedfile_dfs = pl.concat(normalized_bedfiles)
        satellite_db = extract_regions_by_name(
            bedfile_dfs, args.fasta, args.kmer, args.verbose
        )
        # histogram_db = extract_histograms_by_name(bedfile_dfs,args.fasta,args.kmer,True)

        if satellite_metadata:
            for key in satellite_metadata:
                subtypes = [s.lower() for s in satellite_metadata[key]]
                print(subtypes)

                # Filter subtypes that exist in the database
                sub_kmer_db = {
                    name: satellite_db[name.lower()]
                    for name in subtypes
                    if name.lower() in satellite_db
                }

                if not sub_kmer_db:
                    continue  # Skip if none of the subtypes are found

                outputprefix = f"{key}.db"
                outputname = os.path.join(args.directory, outputprefix)
                for subtype_key, value in sub_kmer_db.items():
                    print(
                        f"Saving {subtype_key} into {outputname} containing {len(value)} k-mers"
                    )
                os.makedirs(args.directory, exist_ok=True)
                save_kmer_sets_shared_k(
                    sub_kmer_db, k=args.kmer, output_path=outputname
                )

        else:
            for key, value in satellite_db.items():
                subtypes = [key]
                sub_kmer_db = {name: satellite_db[name] for name in subtypes}
                outputprefix = f"{key}.db"
                outputname = os.path.join(args.directory, outputprefix)
                print(f"Saving {outputname} containing {len(value)} k-mers")
                os.makedirs(args.directory, exist_ok=True)
                save_kmer_sets_shared_k(
                    sub_kmer_db, k=args.kmer, output_path=outputname
                )

    # -------- ANNOTATE LOGIC --------#
    elif args.command == "annotate":
        missing_fastas = [fasta for fasta in args.fasta if not os.path.isfile(fasta)]
        if missing_fastas:
            for fasta in missing_fastas:
                print(f"[ERROR] FASTA file does not exist: {fasta}", file=sys.stderr)
            sys.exit(1)

        # Prep args
        band_height = int(args.band * 1_000_000)
        try:
            windows = derive_window_sizes(args.window)
        except ValueError as error:
            print(f"[ERROR] Invalid --window values: {error}.", file=sys.stderr)
            raise SystemExit(1)
        win = windows[0]
        additional_windows = windows[1:]
        interval = (win + 1) // 2
        directory = args.directory or os.getcwd()
        hash_cache_dir = args.cache_dir

        if not args.quiet:
            label_width = 20

            print(f"{'Output directory:':<{label_width}} {directory}")
            print(f"{'K-mer length:':<{label_width}} {args.kmer}")
            print(f"{'Sketch modulo:':<{label_width}} {args.sketch}")
            print(f"{'Band height:':<{label_width}} {band_height} bp")
            print(f"{'Plot matrices:':<{label_width}} {args.plot}")
            print(f"{'Distal links:':<{label_width}} {args.distal}")
            print(f"{'Threads:':<{label_width}} {args.threads}")
            print(
                f"{'Window sizes:':<{label_width}} "
                f"{', '.join(str(window) for window in windows)} bp"
            )
            print(f"{'ANI threshold:':<{label_width}} {args.identity} %")
            print(
                f"{'Hash cache:':<{label_width}} "
                f"{hash_cache_dir if hash_cache_dir is not None else 'disabled'}"
            )
            print(
                f"{'K-mer dir:':<{label_width}} {args.classify if args.classify else 'None provided'}"
            )
            print("─" * 65)

        # Build list of (fasta, [seq_ids]) pairs
        headers = get_input_headers(args.fasta)
        if args.seq_id:
            pairs = []
            for sid in args.seq_id:
                matches = [f for f, ids in headers if sid in ids]
                if matches:
                    pairs.append((matches[0], [sid]))
                else:
                    # Check if there's sequence bounds for this seq_id
                    seq_bounds = define_bounds(sid)
                    if seq_bounds:
                        print(seq_bounds)
                        matches = [f for f, ids in headers if seq_bounds[0] in ids]
                        if matches:
                            pairs.append((matches[0], [sid]))
                        else:
                            if not args.quiet:
                                print(f"Unable to locate {seq_bounds}. Skipping…")
                    else:
                        if not args.quiet:
                            print(f"Unable to locate {sid}. Skipping…")
        else:
            pairs = headers

        # 3) Open all FASTAs once
        fasta_handles = {f: pysam.FastaFile(f) for f, _ in pairs}
        try:
            satellite_dsu = SatelliteDSU()
            # cache everything into locals
            k_param = args.kmer
            verbosity = args.verbose
            build_sets = build_kmer_sets
            imat = intersection_matrix_thresholded
            get_span = get_diagonal_span
            merge_intv = merge_shared_boundaries

            # 4) Main loops
            for fasta, seq_ids in pairs:
                fh = fasta_handles[fasta]
                for seq_id in seq_ids:
                    sequence_started = time.perf_counter()
                    stage_times = {
                        "hash_wait": 0.0,
                        "sketch": 0.0,
                        "scan": 0.0,
                        "candidates": 0.0,
                        "refinement": 0.0,
                        "output": 0.0,
                    }
                    # Define data structure for satellite coordinates. Variable names for sequence name, length, and if samtools was used for coordinates
                    satellite_coordinate_list = []
                    additional_window_candidates = []
                    periodic_lag_proposals = []
                    distal_link_list = []
                    distal_accumulator = (
                        CandidateNeighborhoodAccumulator(
                            window=win,
                            k=k_param,
                            identity=args.identity,
                            halo_windows=2,
                            skip_adjacent_groups=True,
                        )
                        if args.distal
                        else None
                    )
                    seq_len = fh.get_reference_length(seq_id)
                    seq_bounds = define_bounds(seq_id)
                    cached_sequence_hashes = load_cached_sequence_hashes(
                        hash_cache_dir,
                        fasta,
                        seq_id,
                        seq_len,
                        k_param,
                    )

                    if seq_bounds and not args.quiet:
                        print(f"Found bounds for {seq_id}: {seq_bounds}\n")

                    announce_matrix_creation(
                        seq_id,
                        seq_len,
                        cache_hit_dir=(
                            hash_cache_dir
                            if cached_sequence_hashes is not None
                            else None
                        ),
                    )

                    if seq_len < k_param:
                        if not args.quiet:
                            print(
                                f"Sequence {seq_id} is shorter than k={k_param}; skipping.\n"
                            )
                        continue

                    # If seq_len is close to the band height, then keep everything in the same band.
                    sequence_band_height = band_height
                    if seq_len < band_height + (band_height // 2):
                        sequence_band_height = seq_len
                        if not args.quiet:
                            print(f"Adjusting band to {sequence_band_height} bp.\n")

                    band_plan = SequenceBandPlan(
                        sequence_length=seq_len,
                        kmer=k_param,
                        band_height=sequence_band_height,
                        windows=windows,
                    )
                    n_windows = band_plan.band_count
                    max_len = band_plan.window_plans[win].max_len
                    matrix_threads, use_hash_process = annotation_thread_allocation(
                        args.threads,
                        cache_hit=cached_sequence_hashes is not None,
                        band_count=band_plan.band_count,
                    )
                    set_num_threads(matrix_threads)
                    band_iterator = iter(
                        iter_hashed_fasta_bands(
                            fasta,
                            seq_id,
                            band_plan,
                            use_process=use_hash_process,
                            cache_dir=hash_cache_dir,
                            worker_threads=1,
                        )
                    )
                    try:
                        hash_wait_started = time.perf_counter()
                        first_band = next(band_iterator)
                        stage_times["hash_wait"] += (
                            time.perf_counter() - hash_wait_started
                        )
                    except StopIteration:
                        stage_times["hash_wait"] += (
                            time.perf_counter() - hash_wait_started
                        )
                        set_num_threads(args.threads)
                        continue

                    progress_context = alive_bar(
                        n_windows,
                        title=f"Annotating {seq_id}",
                        unit=" matrix",
                        disable=args.quiet,
                        file=sys.stdout,
                    )
                    matrix_progress = progress_context.__enter__()

                    # Create initial window
                    kmers_list = band_plan.hashes_for_window(first_band, win)
                    first_band_sketches = None
                    if additional_windows:
                        first_band_sketches, sketch_runtime = (
                            build_band_window_sketches(
                                first_band, band_plan, windows, sketch=args.sketch
                            )
                        )
                        prev_ov, prev_nov = first_band_sketches[win]
                        stage_times["sketch"] += sketch_runtime
                    else:
                        sketch_started = time.perf_counter()
                        prev_ov, prev_nov = build_sets(
                            kmers_list, max_len, win, interval, sketch=args.sketch
                        )
                        stage_times["sketch"] += (
                            time.perf_counter() - sketch_started
                        )
                    matrix_started = time.perf_counter()
                    if args.plot:
                        (
                            initial_identity_matrix,
                            scan_threshold_matrix,
                        ) = intersection_matrix_with_threshold(
                            prev_ov,
                            prev_nov,
                            k_param,
                            args.identity,
                        )
                        initial_matrix = initial_identity_matrix >= args.identity
                    elif args.distal:
                        initial_identity_matrix = None
                        initial_matrix = imat(
                            prev_ov, prev_nov, k_param, args.identity
                        )
                        scan_threshold_matrix = initial_matrix
                    else:
                        initial_identity_matrix = None
                        initial_matrix = None
                        scan_threshold_matrix = None
                    spans, periodic_proposals = scan_primary_window(
                        prev_ov,
                        prev_nov,
                        win,
                        k_param,
                        args.identity,
                        threshold_matrix=scan_threshold_matrix,
                    )
                    for proposal in periodic_proposals:
                        periodic_lag_proposals.append(
                            PeriodicLagCandidate(
                                start=proposal.start + first_band.start,
                                end=proposal.end + first_band.start,
                                count=proposal.count,
                                period_bp=proposal.period_bp,
                                harmonic_lags=proposal.harmonic_lags,
                                coverage=proposal.coverage,
                                contrast=proposal.contrast,
                            )
                        )
                    stage_times["scan"] += time.perf_counter() - matrix_started

                    (
                        first_additional_candidates,
                        additional_sketch_runtime,
                        additional_scan_runtime,
                    ) = scan_additional_window_candidates(
                        first_band,
                        band_plan,
                        additional_windows,
                        kmer=k_param,
                        identity=args.identity,
                        sketch=args.sketch,
                        prebuilt_sketches=first_band_sketches,
                    )
                    additional_window_candidates.extend(
                        first_additional_candidates
                    )
                    stage_times["sketch"] += additional_sketch_runtime
                    stage_times["scan"] += additional_scan_runtime

                    if n_windows > 1:
                        matrix_progress()
                        # TODO: Remove low count spans
                        """for element in spans:
                            print(element, element[1], element[1]*win, element[0][1]-element[0][0])"""

                        # Replace 0 here with start prefix
                        candidates_started = time.perf_counter()
                        initial_candidates = merge_intv(
                            intervals=spans, prefix=0, window=win, verbose=False
                        )
                        satellite_coordinate_list.extend(initial_candidates)
                        previous_candidates = initial_candidates
                        previous_prefix = 0
                        distal_candidates = combine_band_candidates(
                            initial_candidates,
                            first_additional_candidates,
                            first_band.start,
                        )
                        if distal_accumulator is not None:
                            distal_accumulator.add_band(
                                prev_ov,
                                prev_nov,
                                distal_candidates,
                                prefix=0,
                                threshold_matrix=initial_matrix,
                            )
                            distal_link_list.extend(
                                detect_candidate_to_all_links_from_full_matrix(
                                    initial_matrix,
                                    distal_candidates,
                                    0,
                                    win,
                                )
                            )
                        plot_links = []
                        if args.distal:
                            plot_links = detect_matrix_distal_links(
                                initial_matrix,
                                first_band.start,
                                win,
                                distal_candidates,
                            )
                            distal_link_list.extend(plot_links)
                        if args.plot:
                            heatmap_path = save_matrix_heatmap(
                                initial_matrix,
                                directory,
                                seq_id,
                                1,
                                first_band.start,
                                win,
                                args.identity,
                                distal_links=plot_links,
                                satellite_ranges=candidate_genomic_ranges(
                                    distal_candidates, first_band.start
                                ),
                                identity_matrix=initial_identity_matrix,
                            )
                            if verbosity:
                                print(f"Saved matrix heatmap to {heatmap_path}")
                        stage_times["candidates"] += (
                            time.perf_counter() - candidates_started
                        )
                        # Iterate through independently fetched/hash-overlapped bands.
                        # The producer queues one future band while this process
                        # performs the current band's matrix work.
                        for hashed_band in iter_with_stage_timing(
                            band_iterator, stage_times, "hash_wait"
                        ):
                            w = hashed_band.index + 1
                            kmers_list = band_plan.hashes_for_window(hashed_band, win)

                            band_sketches = None
                            if additional_windows:
                                band_sketches, sketch_runtime = (
                                    build_band_window_sketches(
                                        hashed_band,
                                        band_plan,
                                        windows,
                                        sketch=args.sketch,
                                    )
                                )
                                ov, nov = band_sketches[win]
                                stage_times["sketch"] += sketch_runtime
                            else:
                                sketch_started = time.perf_counter()
                                ov, nov = build_sets(
                                    kmers_list,
                                    max_len,
                                    win,
                                    interval,
                                    sketch=args.sketch,
                                )
                                stage_times["sketch"] += (
                                    time.perf_counter() - sketch_started
                                )

                            matrix_started = time.perf_counter()
                            if args.plot:
                                (
                                    updated_identity_matrix,
                                    scan_threshold_matrix,
                                ) = intersection_matrix_with_threshold(
                                    ov,
                                    nov,
                                    k_param,
                                    args.identity,
                                )
                                updated_matrix = (
                                    updated_identity_matrix >= args.identity
                                )
                            elif args.distal:
                                updated_identity_matrix = None
                                updated_matrix = imat(
                                    ov, nov, k_param, args.identity
                                )
                                scan_threshold_matrix = updated_matrix
                            else:
                                updated_identity_matrix = None
                                updated_matrix = None
                                scan_threshold_matrix = None
                            new_spans, periodic_proposals = scan_primary_window(
                                ov,
                                nov,
                                win,
                                k_param,
                                args.identity,
                                threshold_matrix=scan_threshold_matrix,
                            )
                            for proposal in periodic_proposals:
                                periodic_lag_proposals.append(
                                    PeriodicLagCandidate(
                                        start=proposal.start + hashed_band.start,
                                        end=proposal.end + hashed_band.start,
                                        count=proposal.count,
                                        period_bp=proposal.period_bp,
                                        harmonic_lags=proposal.harmonic_lags,
                                        coverage=proposal.coverage,
                                        contrast=proposal.contrast,
                                    )
                                )
                            matrix_runtime = time.perf_counter() - matrix_started
                            stage_times["scan"] += matrix_runtime

                            prefix_amount = hashed_band.start

                            # print(new_spans)
                            """if verbosity:
                                print(f"Current prefix: {prefix_amount}\n")"""
                            candidates_started = time.perf_counter()
                            local_candidates = merge_intv(
                                intervals=new_spans,
                                prefix=0,
                                window=win,
                                verbose=False,
                            )
                            satellite_coordinate_list.extend(
                                (x + prefix_amount, y + prefix_amount, count)
                                for x, y, count in local_candidates
                            )
                            (
                                band_additional_candidates,
                                additional_sketch_runtime,
                                additional_scan_runtime,
                            ) = scan_additional_window_candidates(
                                hashed_band,
                                band_plan,
                                additional_windows,
                                kmer=k_param,
                                identity=args.identity,
                                sketch=args.sketch,
                                prebuilt_sketches=band_sketches,
                            )
                            additional_window_candidates.extend(
                                band_additional_candidates
                            )
                            stage_times["sketch"] += additional_sketch_runtime
                            stage_times["scan"] += additional_scan_runtime
                            if args.distal:
                                bridge = detect_adjacent_band_bridge(
                                    prev_ov,
                                    prev_nov,
                                    previous_candidates,
                                    previous_prefix,
                                    ov,
                                    nov,
                                    local_candidates,
                                    prefix_amount,
                                    win,
                                    k_param,
                                    args.identity,
                                )
                                distal_link_list.extend(bridge.links)
                                incorporate_seam_candidates(
                                    satellite_coordinate_list,
                                    bridge.seam_candidates,
                                    win,
                                )
                                if args.plot and (
                                    bridge.links or bridge.seam_candidates
                                ):
                                    pair_path = save_adjacent_matrix_heatmap(
                                        initial_matrix,
                                        updated_matrix,
                                        bridge.cross_matrix,
                                        directory,
                                        seq_id,
                                        w - 1,
                                        w,
                                        previous_prefix,
                                        win,
                                        args.identity,
                                        distal_links=bridge.links,
                                        seam_candidates=bridge.seam_candidates,
                                        satellite_ranges=(
                                            candidate_genomic_ranges(
                                                previous_candidates,
                                                previous_prefix,
                                            )
                                            + candidate_genomic_ranges(
                                                local_candidates,
                                                prefix_amount,
                                            )
                                        ),
                                        previous_identity_matrix=(
                                            initial_identity_matrix
                                        ),
                                        current_identity_matrix=(
                                            updated_identity_matrix
                                        ),
                                        cross_identity_matrix=intersection_matrix_rectangular(
                                            prev_ov,
                                            prev_nov,
                                            ov,
                                            nov,
                                            k_param,
                                        ),
                                    )
                                    if verbosity:
                                        print(
                                            f"Saved adjacent matrix heatmap to {pair_path}"
                                        )
                            if distal_accumulator is not None or args.plot:
                                distal_candidates = combine_band_candidates(
                                    local_candidates,
                                    band_additional_candidates,
                                    prefix_amount,
                                )
                            else:
                                distal_candidates = local_candidates
                            if distal_accumulator is not None:
                                distal_accumulator.add_band(
                                    ov,
                                    nov,
                                    distal_candidates,
                                    prefix=prefix_amount,
                                    threshold_matrix=updated_matrix,
                                )
                                distal_link_list.extend(
                                    detect_candidate_to_all_links_from_full_matrix(
                                        updated_matrix,
                                        distal_candidates,
                                        prefix_amount,
                                        win,
                                    )
                                )
                            plot_links = []
                            if args.distal:
                                plot_links = detect_matrix_distal_links(
                                    updated_matrix,
                                    prefix_amount,
                                    win,
                                    distal_candidates,
                                )
                                distal_link_list.extend(plot_links)
                            if args.plot:
                                heatmap_path = save_matrix_heatmap(
                                    updated_matrix,
                                    directory,
                                    seq_id,
                                    w,
                                    prefix_amount,
                                    win,
                                    args.identity,
                                    distal_links=plot_links,
                                    satellite_ranges=candidate_genomic_ranges(
                                        distal_candidates, prefix_amount
                                    ),
                                    identity_matrix=updated_identity_matrix,
                                )
                                if verbosity:
                                    print(f"Saved matrix heatmap to {heatmap_path}")
                            stage_times["candidates"] += (
                                time.perf_counter() - candidates_started
                            )
                            # Roll matrices forward
                            initial_matrix = updated_matrix
                            initial_identity_matrix = updated_identity_matrix
                            prev_ov, prev_nov = ov, nov
                            previous_candidates = local_candidates
                            previous_prefix = prefix_amount

                            matrix_progress()

                    else:
                        # No progress bar in this case
                        candidates_started = time.perf_counter()
                        initial_candidates = merge_intv(
                            intervals=spans, prefix=0, window=win, verbose=False
                        )
                        satellite_coordinate_list.extend(initial_candidates)
                        distal_candidates = combine_band_candidates(
                            initial_candidates,
                            first_additional_candidates,
                            first_band.start,
                        )
                        if distal_accumulator is not None:
                            distal_accumulator.add_band(
                                prev_ov,
                                prev_nov,
                                distal_candidates,
                                prefix=0,
                                threshold_matrix=initial_matrix,
                            )
                            distal_link_list.extend(
                                detect_candidate_to_all_links_from_full_matrix(
                                    initial_matrix,
                                    distal_candidates,
                                    0,
                                    win,
                                )
                            )
                        plot_links = []
                        if args.distal:
                            plot_links = detect_matrix_distal_links(
                                initial_matrix,
                                first_band.start,
                                win,
                                distal_candidates,
                            )
                            distal_link_list.extend(plot_links)
                        if args.plot:
                            heatmap_path = save_matrix_heatmap(
                                initial_matrix,
                                directory,
                                seq_id,
                                1,
                                first_band.start,
                                win,
                                args.identity,
                                distal_links=plot_links,
                                satellite_ranges=candidate_genomic_ranges(
                                    distal_candidates, first_band.start
                                ),
                                identity_matrix=initial_identity_matrix,
                            )
                            if verbosity:
                                print(f"Saved matrix heatmap to {heatmap_path}")
                        stage_times["candidates"] += (
                            time.perf_counter() - candidates_started
                        )
                        matrix_progress()

                    # Exhausting the producer commits a newly generated full-sequence
                    # hash cache. A warm cache is simply exhausted already.
                    for _ in iter_with_stage_timing(
                        band_iterator, stage_times, "hash_wait"
                    ):
                        pass
                    progress_context.__exit__(None, None, None)
                    set_num_threads(args.threads)
                    sequence_hashes = cached_sequence_hashes
                    if sequence_hashes is None:
                        sequence_hashes = load_cached_sequence_hashes(
                            hash_cache_dir,
                            fasta,
                            seq_id,
                            seq_len,
                            k_param,
                        )

                    candidates_started = time.perf_counter()
                    conventional_candidates = [
                        WindowCandidate(start, end, count, win)
                        for start, end, count in satellite_coordinate_list
                    ]
                    periodic_candidates = validate_periodic_lag_candidates(
                        fh,
                        seq_id,
                        periodic_lag_proposals,
                        additional_window_candidates + conventional_candidates,
                        window=win,
                        kmer=k_param,
                        verbose=verbosity,
                    )
                    all_window_candidates = (
                        additional_window_candidates
                        + conventional_candidates
                        + periodic_candidates
                    )
                    selected_window_candidates = select_multi_window_candidates(
                        all_window_candidates
                    )
                    filtered = [
                        (candidate.start, candidate.end, candidate.count)
                        for candidate in selected_window_candidates
                    ]
                    candidate_windows = [
                        candidate.window for candidate in selected_window_candidates
                    ]
                    candidate_matrix_support = [
                        candidate_passes_support(candidate)
                        for candidate in selected_window_candidates
                    ]
                    if verbosity:
                        print(
                            f"Found {len(selected_window_candidates)} potential "
                            "candidates"
                        )
                    satellite_coordinate_list = filtered
                    if len(windows) > 1:
                        os.makedirs(directory, exist_ok=True)
                        audit_candidates = selected_window_candidates
                        if seq_bounds:
                            coordinate_offset = int(seq_bounds[1])
                            audit_candidates = [
                                WindowCandidate(
                                    candidate.start + coordinate_offset,
                                    candidate.end + coordinate_offset,
                                    candidate.count,
                                    candidate.window,
                                    candidate.source,
                                )
                                for candidate in selected_window_candidates
                            ]
                        selection_path = write_window_selection(
                            os.path.join(directory, f"{seq_id}_window_selection.tsv"),
                            audit_candidates,
                        )
                        if verbosity:
                            print(f"Saved multi-window selections to {selection_path}")
                    if distal_accumulator is not None:
                        neighborhood = distal_accumulator.build(
                            satellite_coordinate_list
                        )
                        os.makedirs(directory, exist_ok=True)
                        neighborhood_path = os.path.join(
                            directory, f"{seq_id}_distal_neighborhood.npz"
                        )
                        np.savez_compressed(
                            neighborhood_path,
                            matrix=neighborhood.matrix,
                            genomic_window_starts=neighborhood.genomic_window_starts,
                            candidate_ranges=np.asarray(
                                neighborhood.candidate_ranges, dtype=np.int64
                            ).reshape(-1, 2),
                            window=np.int64(win),
                            halo_windows=np.int64(2),
                            identity=np.int64(args.identity),
                        )
                        distal_link_list.extend(
                            detect_distal_links(
                                neighborhood.matrix,
                                neighborhood.genomic_window_starts,
                                neighborhood.candidate_ranges,
                                win,
                            )
                        )
                        if verbosity:
                            print(
                                f"Saved {neighborhood.matrix.shape[0]}-window distal "
                                f"neighborhood matrix to {neighborhood_path}"
                            )
                    if args.distal:
                        distal_link_list = filter_candidate_distal_links(
                            distal_link_list,
                            satellite_coordinate_list,
                            proximity=2 * win,
                        )
                        distal_link_list = deduplicate_distal_links(distal_link_list)
                    stage_times["candidates"] += (
                        time.perf_counter() - candidates_started
                    )
                    # print(satellite_coordinate_list)
                    # sys.exit(0)
                    if seq_bounds:
                        # seq_bounds[0] is the name seq_bounds[1] is the start offset, 2 is the end offset
                        merged_coordinates = [
                            (
                                seq_bounds[0],  # chrom
                                start + int(seq_bounds[1]),
                                end + int(seq_bounds[1]) - 1,
                                seq_bounds[0],  # name
                                0,  # score
                                ".",  # strand
                                start + int(seq_bounds[1]),  # thickStart
                                # end + seq_bounds[1],      # thickEnd
                                end + int(seq_bounds[1]) - 1,  # thickEnd
                                "255,0,0",  # itemRgb
                            )
                            for start, end, _ in filtered
                        ]
                    else:
                        merged_coordinates = [
                            (
                                seq_id,  # chrom
                                start,  # start
                                end,  # end
                                seq_id,  # name
                                0,  # score
                                ".",  # strand
                                start,  # thickStart
                                end,  # thickEnd
                                "255,0,0",  # itemRgb
                            )
                            for start, end, _ in filtered
                        ]
                    """if verbosity:
                        print(merged_coordinates)"""

                    df1 = pl.DataFrame(
                        merged_coordinates,
                        schema=[
                            "#chrom",
                            "start",
                            "end",
                            "name",
                            "score",
                            "strand",
                            "thickStart",
                            "thickEnd",
                            "itemRgb",
                        ],
                        orient="row",
                    )

                    # If we are using a subseqeunce of a larger fasta, we need to adjust the coordinates back to the original reference frame before outputting
                    refinement_started = time.perf_counter()
                    if seq_bounds:
                        # fa, seq_id, seq_len, band, offset, window, k, df: pl.DataFrame, classify, verbose, quiet) -> None:
                        tuple_of_lists = report_borders(
                            fa=fh,
                            seq_id=seq_id,
                            seq_len=seq_len,
                            band=args.band,
                            offset=seq_bounds[1],
                            window=win,
                            k=k_param,
                            df=df1,
                            classify=args.classify,
                            verbose=verbosity,
                            quiet=args.quiet,
                            sequence_hashes=sequence_hashes,
                            candidate_windows=candidate_windows,
                            candidate_matrix_support=candidate_matrix_support,
                        )
                        (
                            new_starts,
                            new_ends,
                            new_names,
                            monomer,
                            periodicity,
                            hor,
                        ) = tuple_of_lists
                        # print(tuple_of_lists)
                    else:
                        tuple_of_lists = report_borders(
                            fa=fh,
                            seq_id=seq_id,
                            seq_len=seq_len,
                            band=args.band,
                            offset=1,
                            window=win,
                            k=k_param,
                            df=df1,
                            classify=args.classify,
                            verbose=verbosity,
                            quiet=args.quiet,
                            sequence_hashes=sequence_hashes,
                            candidate_windows=candidate_windows,
                            candidate_matrix_support=candidate_matrix_support,
                        )
                        (
                            new_starts,
                            new_ends,
                            new_names,
                            monomer,
                            periodicity,
                            hor,
                        ) = tuple_of_lists
                        # print(tuple_of_lists)
                    stage_times["refinement"] += (
                        time.perf_counter() - refinement_started
                    )

                    # Replace the columns in df1
                    output_started = time.perf_counter()
                    coordinate_offset = int(seq_bounds[1]) if seq_bounds else 0
                    resolve_overlapping_satellites(
                        new_starts,
                        new_ends,
                        new_names,
                        monomer,
                        periodicity,
                        hor,
                        verbose=verbosity,
                    )
                    refined_candidates = list(zip(new_starts, new_ends))
                    if args.distal:
                        distal_link_list = filter_candidate_distal_links(
                            distal_link_list,
                            refined_candidates,
                            proximity=2 * win,
                        )
                        distal_link_list = snap_distal_links_to_candidates(
                            distal_link_list,
                            refined_candidates,
                        )
                        distal_link_list = promote_unmatched_distal_candidates(
                            distal_link_list,
                            new_starts,
                            new_ends,
                            new_names,
                            monomer,
                            periodicity,
                            hor,
                            fasta_file=fh,
                            seq_id=seq_id,
                            seq_len=seq_len,
                            window=win,
                            k=k_param,
                            sequence_hashes=sequence_hashes,
                            verbose=verbosity,
                        )
                        resolve_overlapping_satellites(
                            new_starts,
                            new_ends,
                            new_names,
                            monomer,
                            periodicity,
                            hor,
                            verbose=verbosity,
                        )
                        refined_candidates = list(zip(new_starts, new_ends))
                        distal_link_list = filter_candidate_distal_links(
                            distal_link_list,
                            refined_candidates,
                            proximity=2 * win,
                        )
                        distal_link_list = snap_distal_links_to_candidates(
                            distal_link_list,
                            refined_candidates,
                        )
                        tuple_of_lists = (
                            new_starts,
                            new_ends,
                            new_names,
                            monomer,
                            periodicity,
                            hor,
                        )
                    for (
                        satellite_start,
                        satellite_end,
                        satellite_name,
                        satellite_monomer,
                        satellite_periodicity,
                        satellite_is_hor,
                    ) in zip(
                        new_starts,
                        new_ends,
                        new_names,
                        monomer,
                        periodicity,
                        hor,
                    ):
                        satellite_dsu.add_satellite(
                            seq_id,
                            satellite_start + coordinate_offset,
                            satellite_end + coordinate_offset,
                            name=satellite_name,
                            monomer=satellite_monomer,
                            periodicity=satellite_periodicity,
                            is_hor=satellite_is_hor,
                        )
                    if args.distal:
                        for link in distal_link_list:
                            satellite_dsu.union_by_coordinates(
                                seq_id,
                                link.start1 + coordinate_offset,
                                link.end1 + coordinate_offset,
                                seq_id,
                                link.start2 + coordinate_offset,
                                link.end2 + coordinate_offset,
                            )
                        distal_links_path = os.path.join(
                            directory, f"{seq_id}_distal_links.bedpe"
                        )
                        os.makedirs(directory, exist_ok=True)
                        write_distal_links(
                            distal_link_list,
                            seq_id,
                            distal_links_path,
                            coordinate_offset=coordinate_offset,
                        )
                        if not args.quiet:
                            print(
                                f"Saved {len(distal_link_list)} distal satellite "
                                f"link(s) to {distal_links_path}"
                            )
                    output_starts = [
                        start + coordinate_offset for start in new_starts
                    ]
                    output_ends = [end + coordinate_offset for end in new_ends]
                    item_rgb = satellite_dsu.item_rgb_for_annotations(
                        seq_id,
                        output_starts,
                        output_ends,
                        monomer,
                        periodicity,
                        hor,
                    )
                    if seq_bounds:
                        df2 = pl.DataFrame(
                            {
                                "#chrom": [seq_id] * len(new_starts),
                                "start": output_starts,
                                "end": output_ends,
                                "name": new_names,
                                "score": [e for e in monomer],
                                "strand": ["."] * len(new_starts),
                                "thickStart": output_starts,
                                "thickEnd": output_ends,
                                "itemRgb": item_rgb,
                            }
                        )
                    else:
                        df2 = pl.DataFrame(
                            {
                                "#chrom": [seq_id] * len(new_starts),
                                "start": new_starts,
                                "end": new_ends,
                                "name": new_names,
                                "score": [e for e in monomer],
                                "strand": ["."] * len(new_starts),
                                "thickStart": new_starts,
                                "thickEnd": new_ends,
                                "itemRgb": item_rgb,
                            }
                        )

                    output_format = args.output_format.lower()
                    annotation_filename = f"{seq_id}.{output_format}"
                    annotation_path = os.path.join(directory, annotation_filename)
                    os.makedirs(directory, exist_ok=True)
                    if output_format == "bed":
                        df2.write_csv(annotation_path, separator="\t")
                    else:
                        converted = convert_dataframe_format(df2, output_format)
                        with open(annotation_path, "w") as annotation_handle:
                            annotation_handle.write(converted)

                    csvfilename = (
                        f"{seq_id}_summary.csv"
                        if output_format == "csv"
                        else f"{seq_id}.csv"
                    )
                    csvfilepath = os.path.join(directory, csvfilename)
                    write_summary_file(tuple_of_lists, csvfilepath)
                    stage_times["output"] += time.perf_counter() - output_started

                    if verbosity:
                        total_runtime = time.perf_counter() - sequence_started
                        print(f"Stage timings for {seq_id}:")
                        print(
                            f"  hash/cache wait:    {stage_times['hash_wait']:.3f} s"
                        )
                        print(f"  sketch construction:{stage_times['sketch']:9.3f} s")
                        print(f"  diagonal scan:      {stage_times['scan']:9.3f} s")
                        print(
                            f"  candidate handling: {stage_times['candidates']:9.3f} s"
                        )
                        print(
                            f"  boundary refinement:{stage_times['refinement']:9.3f} s"
                        )
                        print(f"  output:             {stage_times['output']:9.3f} s")
                        print(f"  total:              {total_runtime:9.3f} s")

                    print(
                        f"Successfully finished annotating {seq_id} to "
                        f"{annotation_path}\n"
                    )

            satellite_dsu_path = os.path.join(directory, "satellite_dsu.tsv")
            satellite_dsu_text_path = os.path.join(directory, "satellite_dsu.txt")
            os.makedirs(directory, exist_ok=True)
            satellite_dsu.write_tsv(satellite_dsu_path)
            satellite_dsu.write_text(satellite_dsu_text_path)
            if not args.quiet:
                print(
                    f"Saved {len(satellite_dsu.satellites)} satellite DSU member(s) "
                    f"to {satellite_dsu_path} and {satellite_dsu_text_path}"
                )

        except Exception as e:
            print(e)

    # -------- NTRPRISM LOGIC --------#
    elif args.command == "ntrprism":
        run_ntrprism_command(args)
        return 0
