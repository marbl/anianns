import numpy as np
import anianns.parse_matrix as parse_matrix

from anianns.ani_matrix import (
    intersection_len,
    intersection_matrix,
    intersection_matrix_inverted,
    intersection_matrix_rectangular,
    intersection_matrix_rectangular_thresholded,
    intersection_matrix_thresholded,
    intersection_matrix_with_threshold,
)
from anianns.parse_matrix import (
    CandidateNeighborhoodAccumulator,
    DistalSatelliteLink,
    build_candidate_neighborhood_matrix,
    detect_distal_links,
    detect_adjacent_band_bridge,
    detect_candidate_to_all_links_from_matrix,
    detect_periodic_lag_candidates_from_matches,
    detect_periodic_lag_candidates_from_matrix,
    filter_candidate_distal_links,
    get_diagonal_span,
    get_diagonal_span_from_sets,
    snap_distal_links_to_candidates,
    sobel_edge_mask,
    sobel_edge_response,
)


def test_intersection_len_counts_sorted_overlaps():
    assert intersection_len.py_func(np.array([1, 2, 4]), np.array([2, 3, 4])) == 2


def test_intersection_matrix_is_symmetric_and_handles_empty_windows():
    overlapping = [
        np.array([1, 2, 3], dtype=np.int32),
        np.array([2, 3, 4], dtype=np.int32),
        np.array([], dtype=np.int32),
    ]
    non_overlapping = [
        np.array([1, 2], dtype=np.int32),
        np.array([2, 4], dtype=np.int32),
        np.array([], dtype=np.int32),
    ]

    matrix = intersection_matrix.py_func(overlapping, non_overlapping, 2)

    assert matrix.shape == (3, 3)
    assert np.allclose(matrix, matrix.T)
    assert matrix[0, 0] == 100.0
    assert matrix[2, 2] == 0.0
    assert matrix[0, 2] == 0.0
    assert matrix[0, 1] > 0.0


def test_rectangular_identity_matrix_returns_exact_scores():
    left_overlapping = [np.array([1, 2, 3], dtype=np.int32)]
    left_non_overlapping = [np.array([1, 2], dtype=np.int32)]
    right_overlapping = [
        np.array([1, 2], dtype=np.int32),
        np.array([8, 9], dtype=np.int32),
    ]
    right_non_overlapping = [
        np.array([1, 2, 3], dtype=np.int32),
        np.array([8, 9], dtype=np.int32),
    ]

    matrix = intersection_matrix_rectangular.py_func(
        left_overlapping,
        left_non_overlapping,
        right_overlapping,
        right_non_overlapping,
        21,
    )

    assert matrix.shape == (1, 2)
    assert matrix[0, 0] == 100.0
    assert matrix[0, 1] == 0.0


def test_sobel_edge_mask_outlines_matrix_blocks():
    matrix = np.zeros((9, 9), dtype=np.bool_)
    matrix[2:7, 3:6] = True

    edges = sobel_edge_mask(matrix)

    assert edges.dtype == np.bool_
    assert edges.shape == matrix.shape
    assert edges[2, 3]
    assert edges[6, 5]
    assert not edges[4, 4]
    assert not np.any(sobel_edge_mask(np.zeros((3, 3), dtype=np.bool_)))


def test_sobel_edge_response_preserves_normalized_gradient_strength():
    matrix = np.zeros((9, 9), dtype=np.bool_)
    matrix[2:7, 3:6] = True

    response = sobel_edge_response(matrix)

    assert response.dtype == np.float32
    assert response.shape == matrix.shape
    assert np.isclose(response.max(), 1.0)
    assert 0 < response[2, 3] <= 1
    assert response[4, 4] == 0


def test_periodic_lag_detector_recovers_parallel_diagonal_stripes():
    matches = np.zeros((51, 200), dtype=np.bool_)
    for lag in (10, 20, 30, 40, 50):
        matches[lag, 20 : 180 - lag] = True

    [candidate] = detect_periodic_lag_candidates_from_matches(matches, 1000)

    assert candidate.period_bp == 10_000
    assert candidate.start <= 20_000
    assert candidate.end >= 180_000
    assert candidate.coverage >= 0.60
    assert len(candidate.harmonic_lags) >= 3


def test_periodic_lag_detector_rejects_single_stripe_and_dense_background():
    single = np.zeros((51, 200), dtype=np.bool_)
    single[10, 20:170] = True
    assert detect_periodic_lag_candidates_from_matches(single, 1000) == []

    dense = np.zeros((51, 200), dtype=np.bool_)
    for lag in range(3, 51):
        dense[lag, 20 : 180 - lag] = True
    assert detect_periodic_lag_candidates_from_matches(dense, 1000) == []


def test_periodic_lag_detector_reuses_dense_matrix_exactly():
    matches = np.zeros((51, 200), dtype=np.bool_)
    matrix = np.zeros((200, 200), dtype=np.bool_)
    for lag in (10, 20, 30, 40, 50):
        matches[lag, 20 : 180 - lag] = True
        rows = np.arange(20, 180 - lag)
        matrix[rows, rows + lag] = True
        matrix[rows + lag, rows] = True

    assert detect_periodic_lag_candidates_from_matrix(
        matrix, 1000, max_lag_windows=50
    ) == detect_periodic_lag_candidates_from_matches(matches, 1000)


def test_thresholded_matrix_is_compact_and_matches_score_threshold():
    overlapping = [
        np.array([1, 2, 3, 4], dtype=np.int32),
        np.array([1, 8, 9, 10], dtype=np.int32),
    ]
    non_overlapping = [
        np.array([1, 2, 3, 4], dtype=np.int32),
        np.array([1, 8, 9, 10], dtype=np.int32),
    ]

    scores = intersection_matrix.py_func(overlapping, non_overlapping, 2)
    thresholded = intersection_matrix_thresholded.py_func(
        overlapping, non_overlapping, 2, 75
    )

    assert thresholded.dtype == np.bool_
    assert np.array_equal(thresholded, scores >= 75)
    assert np.array_equal(thresholded, thresholded.T)


def test_combined_identity_matrix_preserves_both_existing_results():
    rng = np.random.default_rng(7)
    non_overlapping = [
        np.sort(rng.choice(500, size=size, replace=False)).astype(np.int32)
        for size in (40, 35, 0, 50, 20)
    ]
    overlapping = [
        np.sort(rng.choice(500, size=size, replace=False)).astype(np.int32)
        for size in (60, 45, 0, 70, 30)
    ]

    expected_identity = intersection_matrix.py_func(
        overlapping, non_overlapping, 21
    )
    expected_threshold = intersection_matrix_thresholded.py_func(
        overlapping, non_overlapping, 21, 86
    )
    identity, threshold = intersection_matrix_with_threshold.py_func(
        overlapping, non_overlapping, 21, 86
    )

    assert np.array_equal(identity, expected_identity)
    assert np.array_equal(threshold, expected_threshold)


def test_rectangular_matrix_compares_selected_windows_between_bands():
    first = np.array([1, 2, 3, 4], dtype=np.int32)
    different = np.array([20, 21, 22, 23], dtype=np.int32)
    matrix = intersection_matrix_rectangular_thresholded.py_func(
        [first, different],
        [first, different],
        np.array([0], dtype=np.int64),
        [different, first],
        [different, first],
        np.array([0, 1], dtype=np.int64),
        2,
        75,
    )

    assert matrix.shape == (1, 2)
    assert matrix.tolist() == [[False, True]]


def test_matrix_free_diagonal_scan_matches_dense_matrix_spans():
    overlapping = [
        np.array([1, 2, 3, 4], dtype=np.int32),
        np.array([1, 2, 3, 4], dtype=np.int32),
        np.array([1, 2, 3, 4], dtype=np.int32),
        np.array([8, 9, 10, 11], dtype=np.int32),
    ]
    non_overlapping = [value.copy() for value in overlapping]
    dense = intersection_matrix_thresholded.py_func(
        overlapping, non_overlapping, 2, 75
    )

    assert get_diagonal_span_from_sets(
        overlapping, non_overlapping, 1000, 2, 75, 2
    ) == get_diagonal_span(dense, 1000, 2)


def test_matrix_free_diagonal_scan_matches_dense_randomized_inputs():
    rng = np.random.default_rng(42)
    non_overlapping = [
        np.sort(rng.choice(200, size=30, replace=False)).astype(np.int32)
        for _ in range(20)
    ]
    overlapping = [value.copy() for value in non_overlapping]
    dense = intersection_matrix_thresholded.py_func(
        overlapping, non_overlapping, 3, 70
    )

    assert get_diagonal_span_from_sets(
        overlapping, non_overlapping, 500, 3, 70, 2
    ) == get_diagonal_span(dense, 500, 2)


def test_candidate_neighborhood_matrix_retains_mapping_and_distal_cells():
    overlapping = [
        np.array([1, 2, 3, 4], dtype=np.int32),
        np.array([1, 2, 3, 4], dtype=np.int32),
        np.array([8, 9, 10, 11], dtype=np.int32),
        np.array([1, 2, 3, 4], dtype=np.int32),
        np.array([1, 2, 3, 4], dtype=np.int32),
    ]
    non_overlapping = [value.copy() for value in overlapping]

    neighborhood = build_candidate_neighborhood_matrix(
        overlapping,
        non_overlapping,
        [(0, 1000, 3), (3000, 4000, 3)],
        window=1000,
        k=2,
        identity=75,
        prefix=10_000,
        halo_windows=0,
    )

    assert neighborhood.matrix.shape == (2, 2)
    assert neighborhood.local_window_indices.tolist() == [0, 3]
    assert neighborhood.genomic_window_starts.tolist() == [10_000, 13_000]
    assert neighborhood.candidate_ranges == ((10_000, 11_000), (13_000, 14_000))
    assert neighborhood.matrix[0, 1]


def test_neighborhood_accumulator_links_candidates_across_bands():
    sketch = np.array([1, 2, 3, 4], dtype=np.int32)
    accumulator = CandidateNeighborhoodAccumulator(
        window=1000, k=2, identity=75, halo_windows=0
    )
    accumulator.add_band([sketch], [sketch], [(0, 1000, 3)], prefix=0)
    accumulator.add_band([sketch], [sketch], [(0, 1000, 3)], prefix=10_000)

    neighborhood = accumulator.build([(0, 1000, 3), (10_000, 11_000, 3)])

    assert neighborhood.genomic_window_starts.tolist() == [0, 10_000]
    assert neighborhood.matrix.shape == (2, 2)
    assert neighborhood.matrix[0, 1]


def test_neighborhood_accumulator_reuses_precomputed_band_cells(monkeypatch):
    first = np.array([1, 2, 3, 4], dtype=np.int32)
    second = np.array([20, 21, 22, 23], dtype=np.int32)
    accumulator = CandidateNeighborhoodAccumulator(
        window=1000, k=2, identity=75, halo_windows=0
    )
    accumulator.add_band(
        [first, second],
        [first, second],
        [(0, 2000, 3)],
        prefix=0,
        threshold_matrix=np.array([[True, False], [False, True]]),
    )
    assert accumulator._band_matrices[0][1].shape == (2, 2)
    accumulator.add_band(
        [first],
        [first],
        [(0, 1000, 3)],
        prefix=10_000,
        threshold_matrix=np.array([[True]]),
    )

    neighborhood = accumulator.build([(0, 2000, 3), (10_000, 11_000, 3)])

    assert neighborhood.matrix[0, 0]
    assert not neighborhood.matrix[0, 1]
    assert neighborhood.matrix[0, 2]
    assert neighborhood.matrix[2, 2]


def test_neighborhood_accumulator_can_skip_already_bridged_adjacent_bands():
    sketch = np.array([1, 2, 3, 4], dtype=np.int32)
    accumulator = CandidateNeighborhoodAccumulator(
        window=1000,
        k=2,
        identity=75,
        halo_windows=0,
        skip_adjacent_groups=True,
    )
    for prefix in (0, 10_000, 20_000):
        accumulator.add_band(
            [sketch],
            [sketch],
            [(0, 1000, 3)],
            prefix=prefix,
            threshold_matrix=np.array([[True]]),
        )

    neighborhood = accumulator.build(
        [(0, 1000, 3), (10_000, 11_000, 3), (20_000, 21_000, 3)]
    )

    assert not neighborhood.matrix[0, 1]
    assert not neighborhood.matrix[1, 2]
    assert neighborhood.matrix[0, 2]


def test_adjacent_bridge_finds_candidate_to_unknown_region():
    matching = np.array([1, 2, 3, 4], dtype=np.int32)
    different = np.array([20, 21, 22, 23], dtype=np.int32)
    previous = [matching.copy() for _ in range(12)] + [
        different.copy() for _ in range(18)
    ]
    current = [different.copy() for _ in range(15)] + [
        matching.copy() for _ in range(12)
    ] + [different.copy() for _ in range(3)]

    bridge = detect_adjacent_band_bridge(
        previous,
        previous,
        [(0, 12_000, 12)],
        0,
        current,
        current,
        [],
        30_000,
        1000,
        2,
        75,
        seam_halo_windows=5,
    )

    assert len(bridge.links) == 1
    assert (bridge.links[0].start1, bridge.links[0].end1) == (0, 12_000)
    assert (bridge.links[0].start2, bridge.links[0].end2) == (
        45_000,
        57_000,
    )
    assert bridge.cross_matrix.shape == (30, 30)


def test_adjacent_bridge_finds_candidate_free_seam_continuation():
    matching = np.array([1, 2, 3, 4], dtype=np.int32)
    different = np.array([20, 21, 22, 23], dtype=np.int32)
    previous = [different.copy() for _ in range(8)] + [
        matching.copy() for _ in range(12)
    ]
    current = [matching.copy() for _ in range(12)] + [
        different.copy() for _ in range(8)
    ]

    bridge = detect_adjacent_band_bridge(
        previous,
        previous,
        [],
        0,
        current,
        current,
        [],
        20_000,
        1000,
        2,
        75,
        seam_halo_windows=12,
    )

    assert bridge.links == ()
    assert bridge.seam_candidates == ((8000, 32_000, 24),)


def test_distal_link_detection_scores_known_candidate_pairs():
    matrix = np.eye(12, dtype=np.bool_)
    matrix[1:5, 7:11] = True
    matrix[7:11, 1:5] = True
    starts = np.arange(12, dtype=np.int64) * 1000

    links = detect_distal_links(
        matrix,
        starts,
        ((1000, 5000), (7000, 11_000)),
        1000,
        min_windows=3,
    )

    assert len(links) == 1
    link = links[0]
    assert (link.start1, link.end1, link.start2, link.end2) == (
        1000,
        5000,
        7000,
        11_000,
    )
    assert link.source == "candidate"
    assert link.density == 1.0
    assert link.row_coverage == link.column_coverage == 1.0
    assert link.hit_count == 16


def test_distal_link_component_fallback_finds_unassigned_rectangle():
    matrix = np.eye(12, dtype=np.bool_)
    matrix[1:4, 7:10] = True
    matrix[7:10, 1:4] = True
    starts = np.arange(12, dtype=np.int64) * 1000

    links = detect_distal_links(matrix, starts, (), 1000, min_windows=3)

    assert len(links) == 1
    link = links[0]
    assert (link.start1, link.end1, link.start2, link.end2) == (
        1000,
        4000,
        7000,
        10_000,
    )
    assert link.source == "component"
    assert link.hit_count == 9


def test_distal_component_fallback_uses_find_objects_once(monkeypatch):
    matrix = np.eye(30, dtype=np.bool_)
    matrix[1:5, 20:24] = True
    matrix[20:24, 1:5] = True
    original = parse_matrix.ndimage.find_objects
    calls = []

    def recording_find_objects(labels, max_label=0):
        calls.append((labels.shape, max_label))
        return original(labels, max_label=max_label)

    monkeypatch.setattr(parse_matrix.ndimage, "find_objects", recording_find_objects)

    links = detect_distal_links(
        matrix,
        np.arange(30, dtype=np.int64) * 1000,
        (),
        1000,
        min_windows=3,
    )

    assert len(links) == 1
    assert len(calls) == 1
    assert calls[0][0] == matrix.shape


def test_find_objects_component_bounds_match_naive_full_label_scans(monkeypatch):
    rng = np.random.default_rng(7)
    matrix = np.eye(80, dtype=np.bool_)
    noise_rows = rng.integers(0, 35, size=120)
    noise_columns = rng.integers(45, 80, size=120)
    matrix[noise_rows, noise_columns] = True
    matrix[noise_columns, noise_rows] = True
    matrix[5:18, 55:70] = True
    matrix[55:70, 5:18] = True
    starts = np.arange(80, dtype=np.int64) * 1000
    real_find_objects = parse_matrix.ndimage.find_objects

    def naive_find_objects(labels, max_label=0):
        objects = []
        for label in range(1, max_label + 1):
            rows, columns = np.nonzero(labels == label)
            if len(rows) == 0:
                objects.append(None)
            else:
                objects.append(
                    (
                        slice(int(rows.min()), int(rows.max()) + 1),
                        slice(int(columns.min()), int(columns.max()) + 1),
                    )
                )
        return objects

    monkeypatch.setattr(parse_matrix.ndimage, "find_objects", naive_find_objects)
    expected = detect_distal_links(matrix, starts, (), 1000, min_windows=3)
    monkeypatch.setattr(parse_matrix.ndimage, "find_objects", real_find_objects)
    actual = detect_distal_links(matrix, starts, (), 1000, min_windows=3)

    assert actual == expected


def test_distal_link_component_fallback_rejects_sparse_noise():
    matrix = np.eye(12, dtype=np.bool_)
    matrix[1, 8] = matrix[8, 1] = True
    matrix[3, 10] = matrix[10, 3] = True

    assert detect_distal_links(
        matrix,
        np.arange(12, dtype=np.int64) * 1000,
        (),
        1000,
        min_windows=3,
    ) == []


def test_candidate_to_all_scan_finds_one_known_one_unknown_block():
    matrix = np.zeros((12, 40), dtype=np.bool_)
    matrix[:, :12] = True
    matrix[:, 20:32] = True
    selected_starts = np.arange(12, dtype=np.int64) * 1000
    all_starts = np.arange(40, dtype=np.int64) * 1000

    links = detect_candidate_to_all_links_from_matrix(
        matrix,
        selected_starts,
        all_starts,
        ((0, 12_000),),
        1000,
    )

    assert len(links) == 1
    link = links[0]
    assert (link.start1, link.end1, link.start2, link.end2) == (
        0,
        12_000,
        20_000,
        32_000,
    )
    assert link.source == "candidate_to_all"
    assert link.density == 1.0


def test_candidate_to_all_scan_accepts_short_supported_anchor():
    matrix = np.zeros((5, 40), dtype=np.bool_)
    matrix[:, :5] = True
    matrix[:, 20:32] = True

    links = detect_candidate_to_all_links_from_matrix(
        matrix,
        np.arange(5, dtype=np.int64) * 2000,
        np.arange(40, dtype=np.int64) * 2000,
        ((0, 9000),),
        2000,
    )

    assert len(links) == 1
    assert (links[0].start1, links[0].end1) == (0, 9000)
    assert (links[0].start2, links[0].end2) == (40_000, 64_000)


def test_fragmented_edges_do_not_change_distal_prediction():
    matrix = np.zeros((5, 40), dtype=np.bool_)
    for column in range(20, 35):
        row = column % 5
        matrix[row, column] = True
        matrix[(row + 1) % 5, column] = True

    assert detect_candidate_to_all_links_from_matrix(
        matrix,
        np.arange(5, dtype=np.int64) * 2000,
        np.arange(40, dtype=np.int64) * 2000,
        ((0, 9000),),
        2000,
    ) == []


def test_known_candidate_pair_accepts_one_short_high_confidence_axis():
    matrix = np.eye(40, dtype=np.bool_)
    matrix[2:7, 20:32] = True
    matrix[20:32, 2:7] = True

    links = detect_distal_links(
        matrix,
        np.arange(40, dtype=np.int64) * 2000,
        ((4000, 14_000), (40_000, 64_000)),
        2000,
    )

    assert len(links) == 1
    assert links[0].source == "candidate"
    assert (links[0].start1, links[0].end1) == (4000, 14_000)
    assert (links[0].start2, links[0].end2) == (40_000, 64_000)


def test_candidate_to_all_scan_rejects_short_target_runs():
    matrix = np.zeros((12, 30), dtype=np.bool_)
    matrix[:, 20:25] = True

    assert detect_candidate_to_all_links_from_matrix(
        matrix,
        np.arange(12, dtype=np.int64) * 1000,
        np.arange(30, dtype=np.int64) * 1000,
        ((0, 12_000),),
        1000,
    ) == []


def test_candidate_link_filter_removes_links_from_rejected_arrays():
    matrix = np.eye(14, dtype=np.bool_)
    matrix[1:5, 8:12] = True
    matrix[8:12, 1:5] = True
    [link] = detect_distal_links(
        matrix,
        np.arange(14, dtype=np.int64) * 1000,
        ((1000, 5000), (8000, 12_000)),
        1000,
        min_windows=3,
    )

    assert filter_candidate_distal_links(
        [link], [(1000, 5000, 10)]
    ) == []
    assert filter_candidate_distal_links(
        [link], [(1000, 5000, 10), (8000, 12_000, 10)]
    ) == [link]

    self_component = DistalSatelliteLink(
        1000, 4000, 5000, 7000, "component", 0.5, 1.0, 1.0, 8
    )
    assert filter_candidate_distal_links(
        [self_component], [(1000, 5000, 10)], proximity=1000
    ) == []

    snapped = snap_distal_links_to_candidates(
        [link], [(900, 5100), (7900, 12_100)]
    )
    assert (snapped[0].start1, snapped[0].end1) == (900, 5100)
    assert (snapped[0].start2, snapped[0].end2) == (7900, 12_100)


def test_intersection_matrix_inverted_merges_two_blocks():
    a = np.array([[100.0, 80.0], [80.0, 100.0]])
    b = np.array([[100.0, 75.0], [75.0, 100.0]])
    overlapping_a = [np.array([1, 2], dtype=np.int32), np.array([2, 3], dtype=np.int32)]
    non_overlapping_a = [np.array([1, 2], dtype=np.int32), np.array([2, 3], dtype=np.int32)]
    overlapping_b = [np.array([2, 4], dtype=np.int32), np.array([3, 4], dtype=np.int32)]
    non_overlapping_b = [np.array([2, 4], dtype=np.int32), np.array([3, 4], dtype=np.int32)]

    merged = intersection_matrix_inverted.py_func(
        a,
        b,
        overlapping_a,
        non_overlapping_a,
        overlapping_b,
        non_overlapping_b,
        2,
    )

    assert merged.shape == (4, 4)
    assert np.allclose(np.diag(merged), 100.0)
    assert np.allclose(merged, merged.T)
    assert merged[2, 0] > 0.0
