import polars as pl
import numpy as np

from anianns import refine_boundaries as refine


PRISM_RESULT = (171, False, None, [])


def test_recurrent_kmer_set_requires_three_occurrences():
    hashes = np.repeat(np.arange(30, dtype=np.int32), 3)
    hashes = np.concatenate((hashes, np.array([100, 100, 0], dtype=np.int32)))

    result = refine.recurrent_kmer_set(hashes)

    assert result == set(range(1, 30))
    assert 100 not in result
    assert 0 not in result


def test_find_last_matching_index_clamps_out_of_range_search():
    assert refine.find_last_matching_index([10, 20, 30], {20}, 0, 100) == 1
    assert refine.find_last_matching_index([10, 20, 30], {10}, -20, 1) == 0


def test_fixed_window_search_returns_none_when_no_chunk_meets_threshold():
    border_kmers = list(range(400))
    sparse_matches = set(range(10))

    assert (
        refine.find_target_fixed_window_left(
            sparse_matches, border_kmers, 400, 0, False, step_size=4
        )
        is None
    )
    assert (
        refine.find_target_fixed_window_right(
            sparse_matches, border_kmers, 400, 0, False, step_size=4
        )
        is None
    )


def test_fixed_window_search_counts_repeated_matching_kmers():
    border_kmers = ([0] * 100) + ([7] * 100) + ([0] * 100)

    assert (
        refine.find_target_fixed_window_left(
            {7}, border_kmers, 300, 0, False, step_size=3
        )
        == 1
    )
    assert (
        refine.find_target_fixed_window_right(
            {7}, border_kmers, 300, 0, False, step_size=3
        )
        == 1
    )


def test_boundary_histogram_is_suppressed_when_search_requires_extension(capsys):
    border_kmers = ([7] * 10) + ([0] * 90)

    result = refine.find_target_fixed_window_right(
        {7}, border_kmers, 100, 1000, True, step_size=10
    )

    assert result == 0
    assert capsys.readouterr().out == ""


def test_resolved_boundary_histogram_omits_stop_index(capsys):
    border_kmers = ([0] * 50) + ([7] * 10) + ([0] * 40)

    result = refine.find_target_fixed_window_right(
        {7}, border_kmers, 100, 1000, True, step_size=10
    )

    output = capsys.readouterr().out
    assert result == 5
    assert "#" in output
    assert "Boundary stop index" not in output
    assert "Find right boundary point" not in output


def test_boundary_histogram_lines_never_exceed_80_characters(capsys):
    border_kmers = ([0] * 50) + ([7] * 100) + ([0] * 850)

    result = refine.find_target_fixed_window_right(
        {7},
        border_kmers,
        1000,
        12_345_678,
        True,
        step_size=100,
    )

    lines = capsys.readouterr().out.splitlines()
    assert result == 14
    assert lines
    assert max(map(len, lines)) <= 80


def test_boundary_graph_coordinate_labels_never_exceed_80_characters():
    label = refine.boundary_graph_label_line(
        123_456_789_012_345_678_901,
        987_654_321_098_765_432_109,
        73,
    )

    assert len(label) <= 80


def test_detect_precise_boundaries_preserves_candidate_when_extension_fails(
    monkeypatch,
):
    monkeypatch.setattr(refine, "extract_region", lambda **kwargs: "ACTG" * 100)
    monkeypatch.setattr(refine, "ntr_prism", lambda *args, **kwargs: PRISM_RESULT)
    monkeypatch.setattr(refine, "detect_left_boundary", lambda **kwargs: None)
    monkeypatch.setattr(refine, "detect_right_boundary", lambda **kwargs: None)

    result = refine.detect_precise_boundaries(
        fasta_file="input.fa",
        seq_id="chr1",
        seq_len=1000,
        window=20,
        k=4,
        coordinates=(100, 500),
        verbose=False,
        classify=False,
        previous_coordinates=(0, 1),
    )

    assert result == (100, 500, None, PRISM_RESULT)


def test_precise_boundaries_rejects_even_with_distal_metadata(monkeypatch):
    monkeypatch.setattr(refine, "extract_region", lambda **kwargs: "ACTG" * 100)
    monkeypatch.setattr(
        refine, "ntr_prism", lambda *args, **kwargs: (None, False, None)
    )
    monkeypatch.setattr(refine, "detect_left_boundary", lambda **kwargs: 90)
    monkeypatch.setattr(refine, "detect_right_boundary", lambda **kwargs: 510)

    result = refine.detect_precise_boundaries(
        fasta_file="input.fa",
        seq_id="chr1",
        seq_len=1000,
        window=20,
        k=21,
        coordinates=(100, 500),
        verbose=False,
        classify=False,
        previous_coordinates=(0, 1),
    )

    assert result is None


def test_precise_boundaries_still_rejects_without_distal_override(monkeypatch):
    monkeypatch.setattr(refine, "extract_region", lambda **kwargs: "ACTG" * 100)
    monkeypatch.setattr(
        refine, "ntr_prism", lambda *args, **kwargs: (None, False, None)
    )

    result = refine.detect_precise_boundaries(
        fasta_file="input.fa",
        seq_id="chr1",
        seq_len=1000,
        window=20,
        k=21,
        coordinates=(100, 500),
        verbose=False,
        classify=False,
        previous_coordinates=(0, 1),
    )

    assert result is None


def test_precise_boundaries_accepts_k6_when_user_k_rejects(monkeypatch, capsys):
    k6_result = (68, False, None, [])
    observed_kmers = []
    monkeypatch.setattr(refine, "extract_region", lambda **kwargs: "ACTG" * 100)

    def fake_prism(_sequence, _size, kmer, **_kwargs):
        observed_kmers.append(kmer)
        return k6_result if kmer == 6 else (None, False, None)

    monkeypatch.setattr(refine, "ntr_prism", fake_prism)
    monkeypatch.setattr(refine, "detect_left_boundary", lambda **kwargs: 90)
    monkeypatch.setattr(refine, "detect_right_boundary", lambda **kwargs: 510)

    result = refine.detect_precise_boundaries(
        fasta_file="input.fa",
        seq_id="chr1",
        seq_len=1000,
        window=20,
        k=21,
        coordinates=(100, 500),
        verbose=True,
        classify=False,
        previous_coordinates=(0, 1),
        strong_matrix_evidence=True,
    )

    assert observed_kmers == [21, 6]
    assert result == (90, 510, None, k6_result)
    assert capsys.readouterr().out == (
        "Accepted candidate 100-500:\n"
        "    Strong matrix support\n"
        "    NTRPrism k=6: monomer=68 bp\n"
        "    NTRPrism k=21 (user-selected): rejected\n"
        "    Estimated left boundary: 90\n"
        "    Estimated right boundary: 510\n"
    )


def test_precise_boundaries_uses_k6_instead_of_zero_user_k_result(monkeypatch):
    k6_result = (68, False, None, [])
    monkeypatch.setattr(refine, "extract_region", lambda **kwargs: "ACTG" * 100)
    monkeypatch.setattr(
        refine,
        "ntr_prism",
        lambda _sequence, _size, kmer, **_kwargs: (
            k6_result if kmer == 6 else (0, False, None, [])
        ),
    )
    monkeypatch.setattr(refine, "detect_left_boundary", lambda **kwargs: 90)
    monkeypatch.setattr(refine, "detect_right_boundary", lambda **kwargs: 510)

    result = refine.detect_precise_boundaries(
        fasta_file="input.fa",
        seq_id="chr1",
        seq_len=1000,
        window=20,
        k=21,
        coordinates=(100, 500),
        verbose=False,
        classify=False,
        previous_coordinates=(0, 1),
    )

    assert result == (90, 510, None, k6_result)


def test_strong_matrix_evidence_cannot_override_both_ntr_rejections(
    monkeypatch, capsys
):
    observed_kmers = []
    monkeypatch.setattr(refine, "extract_region", lambda **kwargs: "ACTG" * 100)

    def reject(_sequence, _size, kmer, **_kwargs):
        observed_kmers.append(kmer)
        return None, False, None

    monkeypatch.setattr(refine, "ntr_prism", reject)
    monkeypatch.setattr(refine, "detect_left_boundary", lambda **kwargs: 90)
    monkeypatch.setattr(refine, "detect_right_boundary", lambda **kwargs: 510)

    result = refine.detect_precise_boundaries(
        fasta_file="input.fa",
        seq_id="chr1",
        seq_len=1000,
        window=20,
        k=21,
        coordinates=(100, 500),
        verbose=True,
        classify=False,
        previous_coordinates=(0, 1),
        strong_matrix_evidence=True,
    )

    assert observed_kmers == [21, 6]
    assert result is None
    assert capsys.readouterr().out == (
        "Removed candidate 100-500: NTRPrism rejected at k=21 and k=6; "
        "matrix support was strong.\n"
    )


def test_detect_precise_boundaries_accepts_zero_as_valid_left_boundary(monkeypatch):
    monkeypatch.setattr(refine, "extract_region", lambda **kwargs: "ACTG" * 100)
    monkeypatch.setattr(refine, "ntr_prism", lambda *args, **kwargs: PRISM_RESULT)
    monkeypatch.setattr(refine, "detect_left_boundary", lambda **kwargs: 0)
    monkeypatch.setattr(refine, "detect_right_boundary", lambda **kwargs: 450)

    result = refine.detect_precise_boundaries(
        fasta_file="input.fa",
        seq_id="chr1",
        seq_len=1000,
        window=20,
        k=4,
        coordinates=(100, 500),
        verbose=False,
        classify=False,
        previous_coordinates=(0, 1),
    )

    assert result[:2] == (0, 450)


def test_right_boundary_extends_search_instead_of_returning_none(monkeypatch):
    border_hashes = iter([np.array([1]), np.array([2])])
    monkeypatch.setattr(
        refine, "region_canonical_hashes", lambda *args, **kwargs: next(border_hashes)
    )
    targets = iter([2, 1])
    monkeypatch.setattr(
        refine, "find_target_fixed_window_right", lambda **kwargs: next(targets)
    )
    searched = []

    def fake_find(kmers, *args):
        searched.append(kmers)
        return 0

    monkeypatch.setattr(refine, "find_last_matching_index", fake_find)

    result = refine.detect_right_boundary(
        fasta_file="input.fa",
        array_seq="ARRAY",
        array_seq_size=100,
        seq_id="chr1",
        boundary_point=1000,
        limit=None,
        k=4,
        window=100,
        interval=50,
        boundary_chunk_size=100,
        boundary=2000,
        array_kmer_set={9},
    )

    assert result is not None
    assert searched == [[2]]


def test_right_boundary_searches_left_after_empty_initial_range(monkeypatch):
    border_hashes = iter([np.array([1]), np.array([2])])
    monkeypatch.setattr(
        refine, "region_canonical_hashes", lambda *args, **kwargs: next(border_hashes)
    )
    targets = iter([None, 2])
    observed_sets = []

    def fake_target(**kwargs):
        observed_sets.append(kwargs["array_kmer_set"])
        return next(targets)

    monkeypatch.setattr(refine, "find_target_fixed_window_right", fake_target)
    monkeypatch.setattr(refine, "find_last_matching_index", lambda *args: 250)

    result = refine.detect_right_boundary(
        fasta_file="input.fa",
        array_seq="ARRAY",
        array_seq_size=100,
        seq_id="chr1",
        boundary_point=1000,
        limit=None,
        k=4,
        window=100,
        interval=50,
        boundary_chunk_size=100,
        boundary=2000,
        array_kmer_set={9},
        fallback_array_kmer_set={7},
    )

    assert result == 904
    assert observed_sets == [{9}, {7}]


def test_left_boundary_extension_searches_updated_kmers(monkeypatch):
    border_hashes = iter([np.array([1]), np.array([2])])
    monkeypatch.setattr(
        refine, "region_canonical_hashes", lambda *args, **kwargs: next(border_hashes)
    )
    targets = iter([5, 1])
    monkeypatch.setattr(
        refine, "find_target_fixed_window_left", lambda **kwargs: next(targets)
    )
    searched = []

    def fake_find(kmers, *args):
        searched.append(kmers)
        return 0

    monkeypatch.setattr(refine, "find_last_matching_index", fake_find)

    result = refine.detect_left_boundary(
        fasta_file="input.fa",
        array_seq="ARRAY",
        array_seq_size=100,
        seq_id="chr1",
        boundary_point=1000,
        k=4,
        window=100,
        interval=50,
        boundary_chunk_size=100,
        verbosity=False,
        prev_border_coordinate=0,
        boundary=0,
        array_kmer_set={9},
    )

    assert result is not None
    assert searched == [[2]]


def test_report_borders_drops_candidate_when_refinement_raises(monkeypatch):
    monkeypatch.setattr(
        refine,
        "detect_precise_boundaries",
        lambda **kwargs: (_ for _ in ()).throw(RuntimeError("extension failed")),
    )
    df = pl.DataFrame({"start": [100], "end": [500]})

    result = refine.report_borders(
        fa="input.fa",
        seq_id="chr1",
        seq_len=1000,
        band=2.0,
        offset=0,
        window=20,
        k=4,
        df=df,
        classify=False,
        verbose=False,
        quiet=True,
    )

    assert result == ([], [], [], [], [], [])


def test_report_borders_drops_result_with_zero_ntr_monomer(monkeypatch):
    monkeypatch.setattr(
        refine,
        "detect_precise_boundaries",
        lambda **kwargs: (*kwargs["coordinates"], None, (0, False, None, [])),
    )
    df = pl.DataFrame({"start": [100], "end": [500]})

    result = refine.report_borders(
        fa="input.fa",
        seq_id="chr1",
        seq_len=1000,
        band=2.0,
        offset=0,
        window=20,
        k=21,
        df=df,
        classify=False,
        verbose=False,
        quiet=True,
    )

    assert result == ([], [], [], [], [], [])


def test_report_borders_does_not_merge_unrelated_band_edge_candidates(monkeypatch):
    calls = []

    def fake_detect(**kwargs):
        coordinates = kwargs["coordinates"]
        calls.append(coordinates)
        return (*coordinates, None, PRISM_RESULT)

    monkeypatch.setattr(refine, "detect_precise_boundaries", fake_detect)
    df = pl.DataFrame(
        {
            "start": [1_900_000, 2_500_000],
            "end": [2_000_000, 2_600_000],
        }
    )

    refine.report_borders(
        fa="input.fa",
        seq_id="chr1",
        seq_len=3_000_000,
        band=2.0,
        offset=0,
        window=2000,
        k=21,
        df=df,
        classify=False,
        verbose=False,
        quiet=True,
    )

    assert calls == [(1_900_001, 2_000_001), (2_500_001, 2_600_001)]


def test_report_borders_labels_and_separates_verbose_candidates(
    monkeypatch, capsys
):
    monkeypatch.setattr(
        refine,
        "detect_precise_boundaries",
        lambda **kwargs: (*kwargs["coordinates"], None, PRISM_RESULT),
    )
    df = pl.DataFrame(
        {"start": [100, 10_000], "end": [500, 20_000]}
    )

    refine.report_borders(
        fa="input.fa",
        seq_id="chr1",
        seq_len=30_000,
        band=2.0,
        offset=0,
        window=1000,
        k=21,
        df=df,
        classify=False,
        verbose=True,
        quiet=False,
    )

    output = capsys.readouterr().out
    separator = "=" * refine.BOUNDARY_GRAPH_MAX_WIDTH
    assert output.count(separator) == 2
    assert (
        "Analyzing candidate: left boundary=101, right boundary=501" in output
    )
    assert (
        "Analyzing candidate: left boundary=10001, right boundary=20001"
        in output
    )


def test_report_borders_merges_candidates_that_share_the_same_band_edge(monkeypatch):
    calls = []

    def fake_detect(**kwargs):
        coordinates = kwargs["coordinates"]
        calls.append(coordinates)
        return (*coordinates, None, PRISM_RESULT)

    monkeypatch.setattr(refine, "detect_precise_boundaries", fake_detect)
    df = pl.DataFrame(
        {
            "start": [1_900_000, 2_000_000],
            "end": [2_000_000, 2_100_000],
        }
    )

    refine.report_borders(
        fa="input.fa",
        seq_id="chr1",
        seq_len=3_000_000,
        band=2.0,
        offset=0,
        window=2000,
        k=21,
        df=df,
        classify=False,
        verbose=False,
        quiet=True,
    )

    assert calls == [(1_900_001, 2_100_001)]


def test_report_borders_uses_each_candidates_selected_window(monkeypatch):
    calls = []

    def fake_detect(**kwargs):
        calls.append((kwargs["coordinates"], kwargs["window"]))
        return (*kwargs["coordinates"], None, PRISM_RESULT)

    monkeypatch.setattr(refine, "detect_precise_boundaries", fake_detect)
    df = pl.DataFrame({"start": [100, 10_000], "end": [500, 20_000]})

    refine.report_borders(
        fa="input.fa",
        seq_id="chr1",
        seq_len=30_000,
        band=2.0,
        offset=0,
        window=1000,
        candidate_windows=[1000, 4000],
        k=21,
        df=df,
        classify=False,
        verbose=False,
        quiet=True,
    )

    assert calls == [((101, 501), 1000), ((10_001, 20_001), 4000)]


def test_report_borders_propagates_candidate_matrix_support(monkeypatch):
    observed = []

    def fake_detect(**kwargs):
        observed.append(kwargs["strong_matrix_evidence"])
        return (*kwargs["coordinates"], None, PRISM_RESULT)

    monkeypatch.setattr(refine, "detect_precise_boundaries", fake_detect)
    df = pl.DataFrame({"start": [100, 10_000], "end": [500, 20_000]})

    refine.report_borders(
        fa="input.fa",
        seq_id="chr1",
        seq_len=30_000,
        band=2.0,
        offset=0,
        window=1000,
        candidate_windows=[1000, 4000],
        candidate_matrix_support=[True, False],
        k=21,
        df=df,
        classify=False,
        verbose=False,
        quiet=True,
    )

    assert observed == [True, False]


def test_report_borders_uses_finer_window_when_merging_band_pieces(monkeypatch):
    calls = []

    def fake_detect(**kwargs):
        calls.append((kwargs["coordinates"], kwargs["window"]))
        return (*kwargs["coordinates"], None, PRISM_RESULT)

    monkeypatch.setattr(refine, "detect_precise_boundaries", fake_detect)
    df = pl.DataFrame(
        {"start": [1_900_000, 2_000_000], "end": [2_000_000, 2_100_000]}
    )

    refine.report_borders(
        fa="input.fa",
        seq_id="chr1",
        seq_len=3_000_000,
        band=2.0,
        offset=0,
        window=1000,
        candidate_windows=[4000, 1000],
        k=21,
        df=df,
        classify=False,
        verbose=False,
        quiet=True,
    )

    assert calls == [((1_900_001, 2_100_001), 1000)]


def test_report_borders_preserves_intentional_ntr_rejection(monkeypatch):
    monkeypatch.setattr(refine, "detect_precise_boundaries", lambda **kwargs: None)
    df = pl.DataFrame({"start": [100], "end": [500]})

    result = refine.report_borders(
        fa="input.fa",
        seq_id="chr1",
        seq_len=1000,
        band=0.5,
        offset=0,
        window=20,
        k=4,
        df=df,
        classify=False,
        verbose=False,
        quiet=True,
    )

    assert result == ([], [], [], [], [], [])


def test_report_borders_writes_boolean_hor_value(monkeypatch):
    monkeypatch.setattr(
        refine,
        "detect_precise_boundaries",
        lambda **kwargs: (101, 501, None, (171, True, [2], [342])),
    )
    df = pl.DataFrame({"start": [100], "end": [500]})

    result = refine.report_borders(
        fa="input.fa",
        seq_id="chr1",
        seq_len=1000,
        band=0.5,
        offset=0,
        window=20,
        k=4,
        df=df,
        classify=False,
        verbose=False,
        quiet=True,
    )

    assert result[5] == [True]


def test_ntr_prism_ignores_harmonics_without_a_matching_peak(monkeypatch):
    monkeypatch.setattr(
        refine, "forward_kmer_hashes", lambda *args, **kwargs: np.array([1, 2])
    )
    monkeypatch.setattr(refine, "calculate_hash_distances", lambda values: [10])
    monkeypatch.setattr(refine, "top_n_frequent_distances", lambda values, n: [(10, 50)])
    monkeypatch.setattr(refine, "merge_close_values", lambda values, n: [(10, 50)])
    monkeypatch.setattr(
        refine, "hor_test_local_enrichment", lambda *args, **kwargs: (True, [3])
    )

    result = refine.ntr_prism("ACTG" * 30, 100, 21)

    assert result == (10, True, [3], [])


def test_annotation_ntr_prism_uses_k21_by_default(monkeypatch):
    observed = {}

    def capture_hashes(sequence, kmer):
        observed["kmer"] = kmer
        return np.array([1, 2])

    monkeypatch.setattr(refine, "forward_kmer_hashes", capture_hashes)
    monkeypatch.setattr(refine, "calculate_hash_distances", lambda values: [10])
    monkeypatch.setattr(refine, "top_n_frequent_distances", lambda values, n: [(10, 50)])
    monkeypatch.setattr(refine, "merge_close_values", lambda values, n: [(10, 50)])
    monkeypatch.setattr(
        refine, "hor_test_local_enrichment", lambda *args, **kwargs: False
    )

    refine.ntr_prism("ACTG" * 30, 100)

    assert observed["kmer"] == 21


def test_ntr_prism_uses_dominant_peak_when_shortest_is_not_its_harmonic(
    monkeypatch,
):
    """A frequent 2241-bp repeat must not be mislabeled by a stray 2-bp gap."""
    monkeypatch.setattr(
        refine, "forward_kmer_hashes", lambda *args, **kwargs: np.array([1, 2])
    )
    monkeypatch.setattr(refine, "calculate_hash_distances", lambda values: [2, 2241])
    monkeypatch.setattr(
        refine,
        "top_n_frequent_distances",
        lambda values, n: [(2241, 38_564), (2231, 15_717), (2220, 4_387), (2, 3_373)],
    )
    monkeypatch.setattr(refine, "merge_close_values", lambda values, n: values)
    monkeypatch.setattr(
        refine, "hor_test_local_enrichment", lambda *args, **kwargs: False
    )

    result = refine.ntr_prism("ACTG" * 60_000, 240_000, 21)

    assert result == (2241, False, None, [])


def test_ntr_prism_keeps_shortest_peak_when_dominant_peak_is_its_harmonic(
    monkeypatch,
):
    monkeypatch.setattr(
        refine, "forward_kmer_hashes", lambda *args, **kwargs: np.array([1, 2])
    )
    monkeypatch.setattr(refine, "calculate_hash_distances", lambda values: [171, 342])
    monkeypatch.setattr(
        refine,
        "top_n_frequent_distances",
        lambda values, n: [(342, 20_000), (171, 5_000)],
    )
    monkeypatch.setattr(refine, "merge_close_values", lambda values, n: values)
    monkeypatch.setattr(
        refine, "hor_test_local_enrichment", lambda *args, **kwargs: False
    )

    result = refine.ntr_prism("ACTG" * 10_000, 40_000, 21)

    assert result == (171, False, None, [])


def test_region_hashes_slice_existing_sequence_cache(monkeypatch):
    sequence_hashes = np.arange(100, dtype=np.int32)
    monkeypatch.setattr(
        refine,
        "extract_region",
        lambda **kwargs: (_ for _ in ()).throw(AssertionError("FASTA was fetched")),
    )

    result = refine.region_canonical_hashes(
        "input.fa", "chr1", 10, 20, 4, sequence_hashes
    )

    assert np.array_equal(result, sequence_hashes[10:17])


def test_precise_boundary_reuses_one_core_kmer_set(monkeypatch):
    monkeypatch.setattr(refine, "extract_region", lambda **kwargs: "ACTG" * 100)
    monkeypatch.setattr(refine, "ntr_prism", lambda *args, **kwargs: PRISM_RESULT)
    monkeypatch.setattr(
        refine, "canonical_kmer_hashes", lambda *args, **kwargs: np.array([4, 8, 4])
    )
    seen = []

    def fake_left(**kwargs):
        seen.append(kwargs["array_kmer_set"])
        return 100

    def fake_right(**kwargs):
        seen.append(kwargs["array_kmer_set"])
        return 500

    monkeypatch.setattr(refine, "detect_left_boundary", fake_left)
    monkeypatch.setattr(refine, "detect_right_boundary", fake_right)

    refine.detect_precise_boundaries(
        fasta_file="input.fa",
        seq_id="chr1",
        seq_len=1000,
        window=20,
        k=4,
        coordinates=(100, 500),
        verbose=False,
        classify=False,
        previous_coordinates=(0, 1),
    )

    assert seen[0] is seen[1]
    assert seen[0] == {4, 8}
