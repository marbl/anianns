import polars as pl

from anianns import refine_boundaries as refine


PRISM_RESULT = (171, False, None, [])


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
    extracted = iter(["BORDER1", "BORDER2"])
    monkeypatch.setattr(refine, "extract_region", lambda **kwargs: next(extracted))
    monkeypatch.setattr(
        refine,
        "generate_kmers_from_fasta",
        lambda sequence, k, quiet: iter(
            [9] if sequence == "ARRAY" else ([1] if sequence == "BORDER1" else [2])
        ),
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
    )

    assert result is not None
    assert searched == [[2]]


def test_left_boundary_extension_searches_updated_kmers(monkeypatch):
    extracted = iter(["BORDER1", "BORDER2"])
    monkeypatch.setattr(refine, "extract_region", lambda **kwargs: next(extracted))
    monkeypatch.setattr(
        refine,
        "generate_kmers_from_fasta",
        lambda sequence, k, quiet: iter(
            [9] if sequence == "ARRAY" else ([1] if sequence == "BORDER1" else [2])
        ),
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
    )

    assert result is not None
    assert searched == [[2]]


def test_report_borders_keeps_candidate_when_refinement_raises(monkeypatch):
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

    assert result == ([101], [501], [None], [0], [None], [False])


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
        refine, "generate_kmers_from_fasta_forward_only", lambda **kwargs: iter([1, 2])
    )
    monkeypatch.setattr(refine, "calculate_distances", lambda values: [10])
    monkeypatch.setattr(refine, "top_n_frequent_distances", lambda values, n: [(10, 50)])
    monkeypatch.setattr(refine, "merge_close_values", lambda values, n: [(10, 50)])
    monkeypatch.setattr(
        refine, "hor_test_local_enrichment", lambda *args, **kwargs: (True, [3])
    )

    result = refine.ntr_prism("ACTG" * 30, 100, 21)

    assert result == (10, True, [3], [])
