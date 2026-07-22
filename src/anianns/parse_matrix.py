from collections import Counter
from dataclasses import dataclass
import numpy as np
from scipy import ndimage
from statistics import median
import matplotlib.pyplot as plt

from anianns.ani_matrix import (
    diagonal_span_bounds,
    intersection_matrix_cross_groups_thresholded,
    intersection_matrix_rectangular_thresholded,
    intersection_matrix_thresholded,
    intersection_matrix_selected_thresholded,
    intersection_matrix_selected_vs_all_thresholded,
)
from numba.typed import List as NumbaList


@dataclass(frozen=True)
class CandidateNeighborhoodMatrix:
    """Compact similarity matrix and coordinates for candidate neighborhoods."""

    matrix: np.ndarray
    local_window_indices: np.ndarray
    genomic_window_starts: np.ndarray
    candidate_ranges: tuple


@dataclass(frozen=True)
class DistalSatelliteLink:
    """A supported similarity block linking two genomic intervals."""

    start1: int
    end1: int
    start2: int
    end2: int
    source: str
    density: float
    row_coverage: float
    column_coverage: float
    hit_count: int


@dataclass(frozen=True)
class AdjacentBandBridge:
    """Sparse cross-band evidence retained while two neighboring bands overlap."""

    links: tuple
    seam_candidates: tuple
    cross_matrix: np.ndarray


def _block_statistics(binary, row_indices, column_indices):
    if len(row_indices) == 0 or len(column_indices) == 0:
        return 0.0, 0.0, 0.0, 0
    hit_count = 0
    supported_row_count = 0
    supported_columns = np.zeros(len(column_indices), dtype=np.bool_)
    # Bound temporary memory for very large satellite-pair rectangles.
    for start in range(0, len(row_indices), 256):
        row_chunk = row_indices[start : start + 256]
        block = binary[np.ix_(row_chunk, column_indices)]
        hit_count += int(np.count_nonzero(block))
        supported_row_count += int(np.count_nonzero(np.any(block, axis=1)))
        supported_columns |= np.any(block, axis=0)
    if hit_count == 0:
        return 0.0, 0.0, 0.0, 0
    density = hit_count / (len(row_indices) * len(column_indices))
    row_coverage = supported_row_count / len(row_indices)
    column_coverage = np.count_nonzero(supported_columns) / len(column_indices)
    return density, row_coverage, column_coverage, hit_count


def _range_window_indices(genomic_starts, interval, window):
    start, end = interval
    return np.flatnonzero(
        (genomic_starts < end) & ((genomic_starts + window) > start)
    )


def _ranges_overlap(first, second):
    return max(first[0], second[0]) < min(first[1], second[1])


def _link_confidence(link):
    return min(link.row_coverage, link.column_coverage) + link.density


def _reciprocal_overlap(start1, end1, start2, end2):
    overlap = max(0, min(end1, end2) - max(start1, start2))
    if overlap == 0:
        return 0.0
    return min(overlap / (end1 - start1), overlap / (end2 - start2))


def deduplicate_distal_links(links, minimum_overlap=0.5):
    """Collapse symmetric or near-identical links, preferring candidate calls."""
    ordered = sorted(
        links,
        key=lambda link: (
            link.source != "candidate",
            -_link_confidence(link),
            link.start1,
            link.start2,
        ),
    )
    kept = []
    for link in ordered:
        duplicate = any(
            _reciprocal_overlap(
                link.start1, link.end1, prior.start1, prior.end1
            )
            >= minimum_overlap
            and _reciprocal_overlap(
                link.start2, link.end2, prior.start2, prior.end2
            )
            >= minimum_overlap
            for prior in kept
        )
        if not duplicate:
            kept.append(link)
    return sorted(kept, key=lambda link: (link.start1, link.start2))


def filter_candidate_distal_links(
    links, retained_candidates, minimum_overlap=0.5, proximity=0
):
    """Remove rejected candidate links and within-array component artifacts."""
    ranges = [
        (int(start), int(end))
        for start, end, *_rest in retained_candidates
        if int(end) > int(start)
    ]

    def supported(start, end):
        return any(
            _reciprocal_overlap(start, end, candidate_start, candidate_end)
            >= minimum_overlap
            for candidate_start, candidate_end in ranges
        )

    def nearby_candidates(start, end):
        nearby = set()
        for index, (candidate_start, candidate_end) in enumerate(ranges):
            gap = max(candidate_start - end, start - candidate_end, 0)
            if gap <= proximity:
                nearby.add(index)
        return nearby

    filtered = []
    for link in links:
        if link.source == "candidate" and not (
            supported(link.start1, link.end1)
            and supported(link.start2, link.end2)
        ):
            continue
        if link.source == "component" and (
            nearby_candidates(link.start1, link.end1)
            & nearby_candidates(link.start2, link.end2)
        ):
            continue
        filtered.append(link)
    return filtered


def snap_distal_links_to_candidates(links, candidates, minimum_overlap=0.5):
    """Replace coarse window bounds with the best refined candidate bounds."""
    ranges = [
        (int(start), int(end))
        for start, end, *_rest in candidates
        if int(end) > int(start)
    ]

    def snap(start, end):
        matches = [
            (
                _reciprocal_overlap(start, end, candidate_start, candidate_end),
                candidate_start,
                candidate_end,
            )
            for candidate_start, candidate_end in ranges
        ]
        if not matches:
            return start, end
        overlap, candidate_start, candidate_end = max(matches)
        if overlap < minimum_overlap:
            return start, end
        return candidate_start, candidate_end

    snapped = []
    for link in links:
        first = snap(link.start1, link.end1)
        second = snap(link.start2, link.end2)
        snapped.append(
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
    return deduplicate_distal_links(snapped)


def detect_distal_links(
    matrix,
    genomic_window_starts,
    candidate_ranges,
    window,
    *,
    min_windows=10,
    min_candidate_windows=3,
    min_density=0.10,
    min_axis_coverage=0.50,
    component_min_density=0.15,
    component_min_axis_coverage=0.60,
    diagonal_halo_windows=2,
):
    """
    Detect rectangular off-diagonal similarity blocks.

    Known candidate pairs are scored first. Remaining cells are searched with
    connected components so locally coherent links can still be reported when
    one candidate boundary was incomplete or absent.
    """
    binary = np.asarray(matrix, dtype=np.bool_)
    genomic_starts = np.asarray(genomic_window_starts, dtype=np.int64)
    if binary.ndim != 2 or binary.shape[0] != binary.shape[1]:
        raise ValueError("distal-link detection requires a square matrix")
    if binary.shape[0] != len(genomic_starts):
        raise ValueError("matrix and genomic window coordinates must align")
    if window <= 0:
        raise ValueError("window must be positive")
    if min_candidate_windows <= 0:
        raise ValueError("min_candidate_windows must be positive")
    if binary.size == 0:
        return []

    ranges = sorted(
        {
            (int(start), int(end))
            for start, end in candidate_ranges
            if int(end) > int(start)
        }
    )
    range_indices = [
        _range_window_indices(genomic_starts, interval, window)
        for interval in ranges
    ]
    claimed = np.zeros_like(binary)
    links = []

    # A candidate's own square belongs to its diagonal/self signal.
    for indices in range_indices:
        if len(indices):
            claimed[np.ix_(indices, indices)] = True

    for first_index, first_range in enumerate(ranges):
        first_windows = range_indices[first_index]
        if len(first_windows) < min_candidate_windows:
            continue
        for second_index in range(first_index + 1, len(ranges)):
            second_range = ranges[second_index]
            if _ranges_overlap(first_range, second_range):
                continue
            second_windows = range_indices[second_index]
            if len(second_windows) < min_candidate_windows:
                continue
            density, row_coverage, column_coverage, hit_count = _block_statistics(
                binary, first_windows, second_windows
            )
            if (
                density < min_density
                or row_coverage < min_axis_coverage
                or column_coverage < min_axis_coverage
            ):
                continue
            links.append(
                DistalSatelliteLink(
                    start1=first_range[0],
                    end1=first_range[1],
                    start2=second_range[0],
                    end2=second_range[1],
                    source="candidate",
                    density=density,
                    row_coverage=row_coverage,
                    column_coverage=column_coverage,
                    hit_count=hit_count,
                )
            )
            claimed[np.ix_(first_windows, second_windows)] = True
            claimed[np.ix_(second_windows, first_windows)] = True

    residual = binary & ~claimed
    diagonal_distance = diagonal_halo_windows * window
    for row, genomic_start in enumerate(genomic_starts):
        left = np.searchsorted(
            genomic_starts, genomic_start - diagonal_distance, side="left"
        )
        right = np.searchsorted(
            genomic_starts, genomic_start + diagonal_distance, side="right"
        )
        residual[row, left:right] = False
    residual = np.triu(residual, k=1)

    # Close one-window holes without losing components that touch an edge.
    padded = np.pad(residual, 1, mode="constant")
    joined = ndimage.binary_closing(
        padded, structure=np.ones((3, 3), dtype=np.bool_)
    )[1:-1, 1:-1]
    joined |= residual
    labels, component_count = ndimage.label(
        joined, structure=np.ones((3, 3), dtype=np.int8)
    )

    # ``find_objects`` scans the label image once and returns the exact bounding
    # slice for every label.  Comparing ``labels == label`` here used to rescan
    # the complete matrix once per component, which became quadratic in the
    # number of small noise components on chromosome-scale runs.
    component_objects = ndimage.find_objects(labels, max_label=component_count)
    for component_slice in component_objects:
        if component_slice is None:
            continue
        row_slice, column_slice = component_slice
        row_start, row_stop = int(row_slice.start), int(row_slice.stop)
        column_start, column_stop = int(column_slice.start), int(column_slice.stop)
        if (
            row_stop - row_start < min_windows
            or column_stop - column_start < min_windows
        ):
            continue
        original = residual[row_start:row_stop, column_start:column_stop]
        supported_rows = np.flatnonzero(np.any(original, axis=1)) + row_start
        supported_columns = np.flatnonzero(np.any(original, axis=0)) + column_start
        if (
            len(supported_rows) < min_windows
            or len(supported_columns) < min_windows
        ):
            continue
        density, row_coverage, column_coverage, hit_count = _block_statistics(
            residual, supported_rows, supported_columns
        )
        if (
            density < component_min_density
            or row_coverage < component_min_axis_coverage
            or column_coverage < component_min_axis_coverage
        ):
            continue

        first = (
            int(genomic_starts[supported_rows].min()),
            int(genomic_starts[supported_rows].max()) + window,
        )
        second = (
            int(genomic_starts[supported_columns].min()),
            int(genomic_starts[supported_columns].max()) + window,
        )
        if first[0] > second[0]:
            first, second = second, first
            row_coverage, column_coverage = column_coverage, row_coverage
        if _ranges_overlap(first, second):
            continue
        links.append(
            DistalSatelliteLink(
                start1=first[0],
                end1=first[1],
                start2=second[0],
                end2=second[1],
                source="component",
                density=density,
                row_coverage=row_coverage,
                column_coverage=column_coverage,
                hit_count=hit_count,
            )
        )

    return deduplicate_distal_links(links)


def detect_candidate_to_all_links_from_matrix(
    matrix,
    selected_genomic_starts,
    all_genomic_starts,
    candidate_ranges,
    window,
    *,
    min_windows=10,
    min_anchor_windows=3,
    min_density=0.15,
    min_axis_coverage=0.60,
):
    """Find distal blocks with one axis anchored to a known candidate."""
    binary = np.asarray(matrix, dtype=np.bool_)
    selected_starts = np.asarray(selected_genomic_starts, dtype=np.int64)
    all_starts = np.asarray(all_genomic_starts, dtype=np.int64)
    if binary.shape != (len(selected_starts), len(all_starts)):
        raise ValueError("candidate-to-all matrix coordinates do not align")
    if binary.size == 0:
        return []
    if min_anchor_windows <= 0:
        raise ValueError("min_anchor_windows must be positive")

    links = []
    for candidate_start, candidate_end in sorted(set(candidate_ranges)):
        candidate_rows = _range_window_indices(
            selected_starts, (candidate_start, candidate_end), window
        )
        if len(candidate_rows) < min_anchor_windows:
            continue
        block = binary[candidate_rows].copy()
        self_columns = _range_window_indices(
            all_starts, (candidate_start, candidate_end), window
        )
        block[:, self_columns] = False

        column_support = np.mean(block, axis=0) >= min_axis_coverage
        padded = np.pad(column_support, 1, mode="constant")
        supported = ndimage.binary_closing(
            padded, structure=np.ones(3, dtype=np.bool_)
        )[1:-1]
        supported |= column_support
        changes = np.diff(np.pad(supported.astype(np.int8), 1))
        run_starts = np.flatnonzero(changes == 1)
        run_stops = np.flatnonzero(changes == -1)

        for run_start, run_stop in zip(run_starts, run_stops):
            if run_stop - run_start < min_windows:
                continue
            column_indices = np.arange(run_start, run_stop, dtype=np.int64)
            density, row_coverage, column_coverage, hit_count = _block_statistics(
                binary, candidate_rows, column_indices
            )
            if (
                density < min_density
                or row_coverage < min_axis_coverage
                or column_coverage < min_axis_coverage
            ):
                continue
            target = (
                int(all_starts[run_start]),
                int(all_starts[run_stop - 1]) + window,
            )
            candidate = (int(candidate_start), int(candidate_end))
            if _ranges_overlap(candidate, target):
                continue
            if candidate[0] <= target[0]:
                first, second = candidate, target
            else:
                first, second = target, candidate
                row_coverage, column_coverage = column_coverage, row_coverage
            links.append(
                DistalSatelliteLink(
                    start1=first[0],
                    end1=first[1],
                    start2=second[0],
                    end2=second[1],
                    source="candidate_to_all",
                    density=density,
                    row_coverage=row_coverage,
                    column_coverage=column_coverage,
                    hit_count=hit_count,
                )
            )
    return deduplicate_distal_links(links)


def detect_candidate_to_all_links(
    overlapping,
    non_overlapping,
    candidates,
    prefix,
    window,
    k,
    identity,
):
    """Build and scan a bounded candidate-by-all matrix for standard mode."""
    selected_indices, candidate_ranges = candidate_neighborhood_indices(
        candidates,
        window,
        len(overlapping),
        prefix=prefix,
        halo_windows=0,
    )
    if len(selected_indices) == 0:
        return []
    matrix = intersection_matrix_selected_vs_all_thresholded(
        overlapping,
        non_overlapping,
        selected_indices,
        k,
        identity,
    )
    selected_starts = prefix + (selected_indices * window)
    all_starts = prefix + (np.arange(len(overlapping), dtype=np.int64) * window)
    return detect_candidate_to_all_links_from_matrix(
        matrix,
        selected_starts,
        all_starts,
        candidate_ranges,
        window,
    )


def detect_candidate_to_all_links_from_full_matrix(
    matrix, candidates, prefix, window
):
    """Reuse a plot matrix for candidate-to-all scanning without new ANI work."""
    binary = np.asarray(matrix, dtype=np.bool_)
    selected_indices, candidate_ranges = candidate_neighborhood_indices(
        candidates,
        window,
        binary.shape[0],
        prefix=prefix,
        halo_windows=0,
    )
    if len(selected_indices) == 0:
        return []
    selected_starts = prefix + (selected_indices * window)
    all_starts = prefix + (
        np.arange(binary.shape[0], dtype=np.int64) * window
    )
    return detect_candidate_to_all_links_from_matrix(
        binary[selected_indices],
        selected_starts,
        all_starts,
        candidate_ranges,
        window,
    )


def _cross_candidate_to_all_links(
    row_overlapping,
    row_non_overlapping,
    candidates,
    row_prefix,
    column_overlapping,
    column_non_overlapping,
    column_prefix,
    window,
    k,
    identity,
):
    selected_indices, candidate_ranges = candidate_neighborhood_indices(
        candidates,
        window,
        len(row_overlapping),
        prefix=row_prefix,
        halo_windows=0,
    )
    column_indices = np.arange(len(column_overlapping), dtype=np.int64)
    if len(selected_indices) == 0 or len(column_indices) == 0:
        return selected_indices, np.empty(
            (len(selected_indices), len(column_indices)), dtype=np.bool_
        ), []
    matrix = intersection_matrix_rectangular_thresholded(
        row_overlapping,
        row_non_overlapping,
        selected_indices,
        column_overlapping,
        column_non_overlapping,
        column_indices,
        k,
        identity,
    )
    selected_starts = row_prefix + (selected_indices * window)
    column_starts = column_prefix + (column_indices * window)
    links = detect_candidate_to_all_links_from_matrix(
        matrix,
        selected_starts,
        column_starts,
        candidate_ranges,
        window,
    )
    return selected_indices, matrix, links


def _detect_seam_candidates(
    matrix,
    row_starts,
    column_starts,
    boundary,
    window,
    *,
    min_windows=10,
    min_density=0.15,
    min_axis_coverage=0.60,
    seam_gap_windows=2,
):
    """Find coherent repeats continuing through an adjacent-band boundary."""
    binary = np.asarray(matrix, dtype=np.bool_)
    if binary.size == 0:
        return []
    padded = np.pad(binary, 1, mode="constant")
    joined = ndimage.binary_closing(
        padded, structure=np.ones((3, 3), dtype=np.bool_)
    )[1:-1, 1:-1]
    joined |= binary
    labels, component_count = ndimage.label(
        joined, structure=np.ones((3, 3), dtype=np.int8)
    )
    candidates = []
    for component_slice in ndimage.find_objects(
        labels, max_label=component_count
    ):
        if component_slice is None:
            continue
        row_slice, column_slice = component_slice
        if (
            row_slice.stop - row_slice.start < min_windows
            or column_slice.stop - column_slice.start < min_windows
        ):
            continue
        original = binary[row_slice, column_slice]
        supported_rows = np.flatnonzero(np.any(original, axis=1)) + row_slice.start
        supported_columns = (
            np.flatnonzero(np.any(original, axis=0)) + column_slice.start
        )
        if (
            len(supported_rows) < min_windows
            or len(supported_columns) < min_windows
        ):
            continue
        density, row_coverage, column_coverage, _hit_count = _block_statistics(
            binary, supported_rows, supported_columns
        )
        if (
            density < min_density
            or row_coverage < min_axis_coverage
            or column_coverage < min_axis_coverage
        ):
            continue
        row_start = int(row_starts[supported_rows].min())
        row_end = int(row_starts[supported_rows].max()) + window
        column_start = int(column_starts[supported_columns].min())
        column_end = int(column_starts[supported_columns].max()) + window
        maximum_gap = seam_gap_windows * window
        if boundary - row_end > maximum_gap or column_start - boundary > maximum_gap:
            continue
        candidates.append(
            (
                row_start,
                column_end,
                len(supported_rows) + len(supported_columns),
            )
        )
    return candidates


def detect_adjacent_band_bridge(
    previous_overlapping,
    previous_non_overlapping,
    previous_candidates,
    previous_prefix,
    current_overlapping,
    current_non_overlapping,
    current_candidates,
    current_prefix,
    window,
    k,
    identity,
    *,
    seam_halo_windows=50,
):
    """Link neighboring plot bands without materializing a merged full matrix."""
    previous_count = len(previous_overlapping)
    current_count = len(current_overlapping)
    cross_matrix = np.zeros((previous_count, current_count), dtype=np.bool_)
    links = []

    previous_indices, previous_to_current, previous_links = (
        _cross_candidate_to_all_links(
            previous_overlapping,
            previous_non_overlapping,
            previous_candidates,
            previous_prefix,
            current_overlapping,
            current_non_overlapping,
            current_prefix,
            window,
            k,
            identity,
        )
    )
    if len(previous_indices):
        cross_matrix[previous_indices, :] |= previous_to_current
    links.extend(previous_links)

    current_indices, current_to_previous, current_links = (
        _cross_candidate_to_all_links(
            current_overlapping,
            current_non_overlapping,
            current_candidates,
            current_prefix,
            previous_overlapping,
            previous_non_overlapping,
            previous_prefix,
            window,
            k,
            identity,
        )
    )
    if len(current_indices):
        cross_matrix[:, current_indices] |= current_to_previous.T
    links.extend(current_links)

    halo = max(0, int(seam_halo_windows))
    previous_seam_indices = np.arange(
        max(0, previous_count - halo), previous_count, dtype=np.int64
    )
    current_seam_indices = np.arange(
        0, min(halo, current_count), dtype=np.int64
    )
    seam_candidates = []
    if len(previous_seam_indices) and len(current_seam_indices):
        seam_matrix = intersection_matrix_rectangular_thresholded(
            previous_overlapping,
            previous_non_overlapping,
            previous_seam_indices,
            current_overlapping,
            current_non_overlapping,
            current_seam_indices,
            k,
            identity,
        )
        cross_matrix[np.ix_(previous_seam_indices, current_seam_indices)] |= (
            seam_matrix
        )
        previous_starts = previous_prefix + (previous_seam_indices * window)
        current_starts = current_prefix + (current_seam_indices * window)
        seam_candidates = _detect_seam_candidates(
            seam_matrix,
            previous_starts,
            current_starts,
            current_prefix,
            window,
        )

    return AdjacentBandBridge(
        links=tuple(deduplicate_distal_links(links)),
        seam_candidates=tuple(seam_candidates),
        cross_matrix=cross_matrix,
    )


class CandidateNeighborhoodAccumulator:
    """Retain only candidate-adjacent window sketches across sequence bands."""

    def __init__(
        self,
        window,
        k,
        identity,
        halo_windows=2,
        skip_adjacent_groups=False,
    ):
        if halo_windows < 0:
            raise ValueError("halo_windows must be non-negative")
        self.window = window
        self.k = k
        self.identity = identity
        self.halo_windows = halo_windows
        self.skip_adjacent_groups = bool(skip_adjacent_groups)
        self._windows = {}
        self._band_matrices = {}
        self._next_band_id = 0

    def add_band(
        self,
        overlapping,
        non_overlapping,
        candidates,
        prefix,
        threshold_matrix=None,
    ):
        band_id = self._next_band_id
        self._next_band_id += 1
        local_indices, _ = candidate_neighborhood_indices(
            candidates,
            self.window,
            len(overlapping),
            prefix=prefix,
            halo_windows=self.halo_windows,
        )
        if threshold_matrix is not None:
            threshold_matrix = np.asarray(threshold_matrix, dtype=np.bool_)
            expected_shape = (len(overlapping), len(overlapping))
            if threshold_matrix.shape != expected_shape:
                raise ValueError(
                    "precomputed band matrix shape does not match its sketches"
                )
            retained_indices = local_indices.astype(np.int64, copy=True)
            self._band_matrices[band_id] = (
                retained_indices,
                threshold_matrix[np.ix_(retained_indices, retained_indices)].copy(),
            )
        for local_index in local_indices:
            genomic_start = prefix + (int(local_index) * self.window)
            self._windows.setdefault(
                genomic_start,
                (
                    overlapping[local_index],
                    non_overlapping[local_index],
                    band_id,
                    int(local_index),
                ),
            )

    def build(self, retained_candidates):
        candidate_ranges = tuple(
            (int(start), int(end)) for start, end, _count in retained_candidates
        )
        halo = self.halo_windows * self.window
        retained_windows = []
        for genomic_start, sketches in sorted(self._windows.items()):
            if any(
                genomic_start >= start - halo and genomic_start < end + halo
                for start, end in candidate_ranges
            ):
                retained_windows.append((genomic_start, sketches))

        genomic_starts = np.array(
            [start for start, _sketches in retained_windows], dtype=np.int64
        )
        if not retained_windows:
            matrix = np.empty((0, 0), dtype=np.bool_)
        else:
            overlapping = NumbaList()
            non_overlapping = NumbaList()
            group_ids = np.empty(len(retained_windows), dtype=np.int64)
            local_indices = np.empty(len(retained_windows), dtype=np.int64)
            for retained_index, (
                _start,
                (overlap, non_overlap, band_id, local_index),
            ) in enumerate(retained_windows):
                overlapping.append(overlap)
                non_overlapping.append(non_overlap)
                group_ids[retained_index] = band_id
                local_indices[retained_index] = local_index

            retained_groups = np.unique(group_ids)
            if all(group_id in self._band_matrices for group_id in retained_groups):
                matrix = intersection_matrix_cross_groups_thresholded(
                    overlapping,
                    non_overlapping,
                    group_ids,
                    self.k,
                    self.identity,
                    self.skip_adjacent_groups,
                )
                for group_id in retained_groups:
                    positions = np.flatnonzero(group_ids == group_id)
                    source_indices = local_indices[positions]
                    stored_indices, source = self._band_matrices[int(group_id)]
                    source_positions = np.searchsorted(
                        stored_indices, source_indices
                    )
                    if (
                        np.any(source_positions >= len(stored_indices))
                        or not np.array_equal(
                            stored_indices[source_positions], source_indices
                        )
                    ):
                        raise ValueError(
                            "retained candidate window is missing from its band matrix"
                        )
                    matrix[np.ix_(positions, positions)] = source[
                        np.ix_(source_positions, source_positions)
                    ]
            else:
                # Backward-compatible path for callers without a plot matrix.
                matrix = intersection_matrix_thresholded(
                    overlapping,
                    non_overlapping,
                    self.k,
                    self.identity,
                )

        return CandidateNeighborhoodMatrix(
            matrix=matrix,
            local_window_indices=np.arange(len(genomic_starts), dtype=np.int64),
            genomic_window_starts=genomic_starts,
            candidate_ranges=candidate_ranges,
        )


def append_coordinates(satellite_coordinate_list, prefix, matrix, window, off_diagonal):
    M_diag, M_offdiag = split_diagonal_attached(matrix)
    if off_diagonal:
        merged = merge_shared_boundaries(M_offdiag, prefix, False)
    else:
        merged = merge_shared_boundaries(matrix, prefix, False)
        for coordinates in merged:
            satellite_coordinate_list.append(coordinates)
    return satellite_coordinate_list


def check_same_start_end(pairs, s, window, start):
    lo, hi = s - window, s + window
    if start:
        return [(x, y, count) for x, y, count in pairs if lo <= x <= hi]
    else:
        return [(x, y, count) for x, y, count in pairs if lo <= y <= hi]


def check_new_contained(pairs, start, end):
    return [(x, y, count) for x, y, count in pairs if (x > start) and (y < end)]


def check_new_spans(pairs, start, end):
    return [(x, y, count) for x, y, count in pairs if (x < start) and (y > end)]


def find_in_range_x(data, min_val, max_val):
    return [(i, t) for i, t in enumerate(data) if min_val <= t[0] <= max_val]


def find_in_range_y(data, min_val, max_val):
    return [(i, t) for i, t in enumerate(data) if min_val <= t[1] <= max_val]


def find_contained(data, min_val_x, max_val_x, min_val_y, max_val_y):
    return [
        (i, t) for i, t in enumerate(data) if (min_val_x > t[0]) and (max_val_y < t[1])
    ]


def find_spanning(data, min_val_x, max_val_x, min_val_y, max_val_y):
    return [
        (i, t)
        for i, t in enumerate(data)
        if ((min_val_x < t[0]) and (max_val_y > t[1]))
    ]


def get_diagonal_span(matrix, window, zero_tol):
    n = matrix.shape[0]
    lengths = [0] * n
    coords = [(0, 0)] * n

    for i in range(n):
        row = matrix[i]  # local view of row i
        start = end = i
        cnt = 0

        # scan left of diagonal
        zeros = 0
        for j in range(i - 1, -1, -1):
            if row[j] != 0:
                cnt += 1
                start = j
                zeros = 0
            else:
                zeros += 1
                if zeros == zero_tol:
                    break

        # scan right of diagonal
        zeros = 0
        for j in range(i + 1, n):
            if row[j] != 0:
                cnt += 1
                end = j
                zeros = 0
            else:
                zeros += 1
                if zeros == zero_tol:
                    break

        lengths[i] = cnt
        coords[i] = (start * window, (end * window) + window)

    tuple_counts = Counter(coords)
    sorted_items = sorted(
        ((k, v) for k, v in tuple_counts.items() if v >= 3), key=lambda item: item[0][0]
    )
    return sorted_items


def get_diagonal_span_from_sets(
    overlapping, non_overlapping, window, k, identity, zero_tol
):
    """Match ``get_diagonal_span`` without allocating its dense input matrix."""
    starts, ends = diagonal_span_bounds(
        overlapping, non_overlapping, k, identity, zero_tol
    )
    coordinates = [
        (int(start) * window, (int(end) * window) + window)
        for start, end in zip(starts, ends)
    ]
    tuple_counts = Counter(coordinates)
    return sorted(
        ((coordinate, count) for coordinate, count in tuple_counts.items() if count >= 3),
        key=lambda item: item[0][0],
    )


def build_candidate_neighborhood_matrix(
    overlapping,
    non_overlapping,
    candidates,
    window,
    k,
    identity,
    *,
    prefix=0,
    halo_windows=2,
):
    """Build an exact compact matrix over candidate windows and their halo."""
    local_indices, candidate_ranges = candidate_neighborhood_indices(
        candidates,
        window,
        len(overlapping),
        prefix=prefix,
        halo_windows=halo_windows,
    )
    if len(local_indices) == 0:
        matrix = np.empty((0, 0), dtype=np.bool_)
    else:
        matrix = intersection_matrix_selected_thresholded(
            overlapping, non_overlapping, local_indices, k, identity
        )

    return CandidateNeighborhoodMatrix(
        matrix=matrix,
        local_window_indices=local_indices,
        genomic_window_starts=prefix + (local_indices * window),
        candidate_ranges=candidate_ranges,
    )


def candidate_neighborhood_indices(
    candidates, window, window_count, *, prefix=0, halo_windows=2
):
    """Select the union of candidate windows and return their genomic ranges."""
    if halo_windows < 0:
        raise ValueError("halo_windows must be non-negative")

    selected = np.zeros(window_count, dtype=np.bool_)
    candidate_ranges = []
    for start, end, _count in candidates:
        start = max(0, int(start))
        end = max(start, int(end))
        first = max(0, (start // window) - halo_windows)
        stop = min(
            window_count,
            ((end + window - 1) // window) + halo_windows,
        )
        selected[first:stop] = True
        candidate_ranges.append((prefix + start, prefix + end))

    local_indices = np.flatnonzero(selected).astype(np.int64, copy=False)
    return local_indices, tuple(candidate_ranges)


def merge_shared_boundaries(intervals, prefix, window, verbose=True):
    """
    intervals: list of ((start, end), count)
    returns: list of ((start+prefix, end+prefix), total_count)
    """
    out = []  # list of (x, y, count)

    chk_same = check_same_start_end  # expects: (pairs_list, val, window, is_start)
    chk_cont = check_new_contained  # expects: (pairs_list, x, y) checks if new sequence is smaller
    chk_span = (
        check_new_spans  # expects: (pairs_list, x, y) checks if new seq is bigger
    )

    for (x, y), count in intervals:
        # Base case, append to out if empty
        if len(out) == 0:
            out.append((x, y, count))
            continue
        else:
            # Initialize start and end window buffers
            start_range = (x - (window * 2), x + (window * 2))
            end_range = (y - (window * 2), y + (window * 2))
            process_entry(
                out=out,
                x=x,
                y=y,
                count=count,
                start_range=start_range,
                end_range=end_range,
                find_in_range_x=find_in_range_x,
                find_in_range_y=find_in_range_y,
                find_contained=find_contained,
                find_spanning=find_spanning,
                verbose=verbose,
            )
            print
    return [(x + prefix, y + prefix, cnt) for x, y, cnt in out]


def process_entry(
    out,
    x,
    y,
    count,
    start_range,
    end_range,
    find_in_range_x,
    find_in_range_y,
    find_contained,
    find_spanning,
    verbose=False,
):
    """
    Chooses the best prior entry to merge with (if any), then updates 'out'.
    Priority: match on both x & y (same index) > x > y > spanning > contained.
    Searches are done lazily to avoid unnecessary work.
    """

    idx = None
    reason = None

    # Try x match first (most common case) and only then check y for intersection
    c_x = find_in_range_x(data=out, min_val=start_range[0], max_val=start_range[1])
    if c_x:
        # Only compute y if x matched; try to find the same index for both
        c_y = find_in_range_y(data=out, min_val=end_range[0], max_val=end_range[1])
        if c_y:
            ix_x = {i for i, _ in c_x}
            ix_y = {i for i, _ in c_y}
            common = sorted(ix_x & ix_y)
            if common:
                idx = common[0]
                reason = "x & y"
            else:
                idx = c_x[0][0]
                reason = "x"
        else:
            idx = c_x[0][0]
            reason = "x"
    else:
        # No x match; try y
        c_y = find_in_range_y(data=out, min_val=end_range[0], max_val=end_range[1])
        if c_y:
            idx = c_y[0][0]
            reason = "y"
        else:
            # Only now attempt the more general/expensive checks
            c_s = find_spanning(
                data=out, min_val_x=x, max_val_x=x, min_val_y=y, max_val_y=y
            )
            if c_s:
                idx = c_s[0][0]
                reason = "spanning"
            else:
                c_c = find_contained(
                    data=out, min_val_x=x, max_val_x=x, min_val_y=y, max_val_y=y
                )
                if c_c:
                    idx = c_c[0][0]
                    reason = "contained"

    if idx is not None:
        if verbose:
            print(f"Matches previous entry on {reason}")
        _update_out(out, idx, x, y, count, verbose=verbose)
    else:
        out.append((x, y, count))


def sobel_spans(M, prefix):
    spans = []
    i = 0
    while i < len(M):
        # print(M[i])

        h = int(abs(M[i][3]))  # window length from height

        if h == 0 or h > 1900:
            i += 1
            continue

        end = min(i + h, len(M))  # clamp to array length
        if end <= i:  # safety, though h>0 makes this unlikely
            i += 1
            continue

        window = [abs(M[k][3]) for k in range(i, end)]
        med = median(window)

        if 0.75 <= med / h <= 1.25:
            # print("Yes")
            spans.append((i - 1, i + h + 1))
            i += h
        else:
            # print("No")
            i += 1
    spans = [(prefix * a, prefix * b) for (a, b) in spans]
    return spans


def sobel_with_diagonal_probes(M, thresh=0.7, min_thick=1):
    """
    Compute Sobel edges, display the binary edge map, and for each diagonal
    position (i, i) that is empty (False), draw a red vertical line that
    extends up and down until it hits a vertical run of True pixels with
    thickness >= min_thick.

    Parameters
    ----------
    M : 2D array
        Image/matrix to edge-detect.
    thresh : float
        Threshold on normalized Sobel magnitude to make binary_edges.
    min_thick : int
        Minimum contiguous thickness (in pixels) of a vertical edge to stop.
    figsize : tuple
        Matplotlib figure size.
    """
    # --- Sobel edges ---
    sobel_x = ndimage.sobel(M, axis=1)
    sobel_y = ndimage.sobel(M, axis=0)
    sobel_mag = np.hypot(sobel_x, sobel_y)
    max_val = np.max(sobel_mag)
    if max_val > 0:
        sobel_mag = sobel_mag / max_val
    binary_edges = sobel_mag > thresh

    H, W = binary_edges.shape
    N = min(H, W)

    def stop_y(y0, x, dy):
        """
        March from y0 in direction dy (+1 down, -1 up) until:
          - we hit image border, or
          - we encounter a vertical run of True pixels with length >= min_thick
            at column x, starting at the next step in the marching direction.
        Returns the last y BEFORE the blocking run/border.
        """
        y = y0
        while True:
            ny = y + dy
            if ny < 0 or ny >= H:
                return y  # hit border

            # Check if a vertical run with length >= min_thick begins at ny
            if dy > 0:
                end = min(ny + min_thick, H)
                if end - ny == min_thick and np.all(binary_edges[ny:end, x]):
                    return y
            else:  # dy < 0
                start = max(ny - (min_thick - 1), 0)
                if ny - start + 1 >= min_thick and np.all(
                    binary_edges[start : ny + 1, x]
                ):
                    return y

            y = ny  # keep marching

    # --- Plot base image and mask ---
    plt.figure(figsize=(8, 8))
    plt.imshow(binary_edges, cmap="gray_r", interpolation="nearest")
    plt.title("Sobel Edge Magnitude with Diagonal Probes")
    plt.colorbar(label="Edge (binary)")

    # --- For each diagonal position that is empty, drop a probe line ---
    for i in range(N):
        if not binary_edges[i, i]:  # "space on the diagonal"
            y_top = stop_y(i, i, dy=-1)
            y_bot = stop_y(i, i, dy=+1)
            # Draw the vertical line at x=i from y_top to y_bot
            plt.plot([i, i], [y_top, y_bot], "-", linewidth=1.5, color="red", alpha=0.9)

    plt.tight_layout()
    plt.show()
    ranges = []
    for i in range(N):
        if not binary_edges[i, i]:  # "space on the diagonal"
            y_top = stop_y(i, i, dy=-1)
            y_bot = stop_y(i, i, dy=+1)
            ranges.append((i, y_top, y_bot, y_top - y_bot))
    return ranges


def split_diagonal_attached(M):
    binary = M > 0
    vertical_structure = np.array([[0, 1, 0], [0, 1, 0], [0, 1, 0]], dtype=int)
    horizontal_structure = np.array([[0, 0, 0], [1, 1, 1], [0, 0, 0]], dtype=int)
    vertical_labeled, _ = ndimage.label(binary, structure=vertical_structure)
    horizontal_labeled, _ = ndimage.label(binary, structure=horizontal_structure)

    diag_indices = np.arange(min(M.shape))
    vertical_diag_labels = np.unique(vertical_labeled[diag_indices, diag_indices])
    horizontal_diag_labels = np.unique(horizontal_labeled[diag_indices, diag_indices])

    vertical_keep_mask = np.isin(vertical_labeled, vertical_diag_labels)
    horizontal_keep_mask = np.isin(horizontal_labeled, horizontal_diag_labels)

    combined_mask = np.minimum(vertical_keep_mask, horizontal_keep_mask)

    M_diag = M * combined_mask
    M_offdiag = M * (~combined_mask)

    return M_diag, M_offdiag


def _update_out(out, index, x, y, count, verbose=False):
    old_x, old_y, old_count = out[index]
    min_x = min(old_x, x)
    max_y = max(old_y, y)
    # Keep old bounds if the new segment is "small" vs existing
    if count < old_count / 2:
        out[index] = (old_x, old_y, old_count + count)
    else:
        out[index] = (min_x, max_y, old_count + count)
    if verbose:
        print(f"Index: {index}")
        print(f"Updated to {out[index]}\n")
