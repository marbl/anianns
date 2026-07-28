import csv
import sys
from pathlib import Path

import polars as pl
import pytest

from anianns import anianns as cli
from anianns.kmer_pipeline import HashedBand
from anianns.parse_matrix import DistalSatelliteLink


class FakeFastaFile:
    def __init__(self, sequences):
        self._sequences = sequences
        self.references = list(sequences)

    def fetch(self, seq_id):
        return self._sequences[seq_id]

    def get_reference_length(self, seq_id):
        return len(self._sequences[seq_id])

    def close(self):
        return None


def test_mask_type_and_parser_defaults():
    assert cli.mask_type("12") == 12
    assert cli.mask_type("alpha") == "alpha"

    parser = cli.get_parser()
    args = parser.parse_args(["annotate", "-f", "input.fa"])

    assert args.command == "annotate"
    assert args.kmer == 21
    assert args.sketch == 4
    assert args.window == 2000
    assert args.plot is False
    assert args.distal is False
    assert args.threads == cli.MAX_THREADS
    assert not hasattr(args, "distal_halo")
    assert args.output_format == "bed"

    dense_args = parser.parse_args(
        ["annotate", "-f", "input.fa", "--sketch", "2"]
    )
    assert dense_args.sketch == 2
    assert args.band == 2.0

    requested_threads = min(2, cli.MAX_THREADS)
    threaded_args = parser.parse_args(
        ["annotate", "-f", "input.fa", "--threads", str(requested_threads)]
    )
    assert threaded_args.threads == requested_threads

    with pytest.raises(SystemExit):
        parser.parse_args(["annotate", "-f", "input.fa", "--threads", "0"])

    multi_args = parser.parse_args(
        ["annotate", "-f", "input.fa", "-w", "5000"]
    )
    assert cli.derive_window_sizes(multi_args.window) == (2500, 5000, 10_000)

    with pytest.raises(SystemExit):
        parser.parse_args(
            ["annotate", "-f", "input.fa", "-w", "1000", "2000"]
        )


def test_format_matrix_runtime_compares_with_previous_matrix():
    assert cli.format_matrix_runtime(1, 3, 1.0) == (
        "Matrix 1/3 completed in 1.000 s (comparison baseline)."
    )
    assert cli.format_matrix_runtime(2, 3, 1.25, 1.0) == (
        "Matrix 2/3 completed in 1.250 s — 0.250 s slower "
        "(+25.0%) than matrix 1."
    )
    assert cli.format_matrix_runtime(3, 3, 1.0, 1.25) == (
        "Matrix 3/3 completed in 1.000 s — 0.250 s faster "
        "(-20.0%) than matrix 2."
    )


def test_annotation_thread_allocation_preserves_total_budget():
    assert cli.annotation_thread_allocation(8, cache_hit=False, band_count=3) == (
        7,
        True,
    )
    assert cli.annotation_thread_allocation(1, cache_hit=False, band_count=3) == (
        1,
        False,
    )
    assert cli.annotation_thread_allocation(8, cache_hit=True, band_count=3) == (
        8,
        False,
    )
    assert cli.annotation_thread_allocation(8, cache_hit=False, band_count=1) == (
        8,
        False,
    )


def test_distal_candidates_include_supported_calls_from_other_windows():
    combined = cli.combine_band_candidates(
        [(100_000, 140_000, 18)],
        [cli.WindowCandidate(6_550_000, 6_559_000, 7, 1000)],
        prefix=6_000_000,
    )

    assert combined == [
        (100_000, 140_000, 18),
        (550_000, 559_000, 7),
    ]


def test_write_distal_links_emits_bedpe_metrics_and_offset(tmp_path):
    output_path = tmp_path / "links.bedpe"
    cli.write_distal_links(
        [
            DistalSatelliteLink(
                100,
                500,
                1000,
                1800,
                "candidate",
                0.75,
                0.8,
                0.9,
                42,
            )
        ],
        "chr1",
        output_path,
        coordinate_offset=10_000,
    )

    rows = list(csv.reader(output_path.open(), delimiter="\t"))
    assert rows[0][:6] == [
        "#chrom1",
        "start1",
        "end1",
        "chrom2",
        "start2",
        "end2",
    ]
    assert rows[1][:11] == [
        "chr1",
        "10100",
        "10500",
        "chr1",
        "11000",
        "11800",
        "distal_link_0001",
        "750",
        ".",
        ".",
        "candidate",
    ]


def test_promote_unmatched_distal_candidate_after_ntr_validation(monkeypatch):
    starts = [100]
    ends = [500]
    names = ["known"]
    monomers = [171]
    periodicities = [None]
    hor_flags = [False]
    link = DistalSatelliteLink(
        100, 500, 1000, 1400, "candidate_to_all", 0.8, 0.9, 0.9, 100
    )
    monkeypatch.setattr(
        cli,
        "detect_precise_boundaries",
        lambda **kwargs: (950, 1450, None, (10, False, None, [])),
    )

    promoted_links = cli.promote_unmatched_distal_candidates(
        [link],
        starts,
        ends,
        names,
        monomers,
        periodicities,
        hor_flags,
        fasta_file="input.fa",
        seq_id="chr1",
        seq_len=2000,
        window=100,
        k=21,
    )

    assert list(zip(starts, ends)) == [(100, 500), (950, 1450)]
    assert names == ["known", None]
    assert monomers == [171, 10]
    assert len(promoted_links) == 1
    assert (promoted_links[0].start2, promoted_links[0].end2) == (950, 1450)


def test_rejects_unmatched_distal_candidate_when_ntr_fails(monkeypatch):
    starts, ends = [100], [500]
    names, monomers, periodicities, hor_flags = [None], [171], [None], [False]
    link = DistalSatelliteLink(
        100, 500, 1000, 1400, "candidate_to_all", 0.8, 0.9, 0.9, 100
    )
    monkeypatch.setattr(cli, "detect_precise_boundaries", lambda **kwargs: None)

    promoted_links = cli.promote_unmatched_distal_candidates(
        [link],
        starts,
        ends,
        names,
        monomers,
        periodicities,
        hor_flags,
        fasta_file="input.fa",
        seq_id="chr1",
        seq_len=2000,
        window=100,
        k=21,
    )

    assert promoted_links == []
    assert list(zip(starts, ends)) == [(100, 500)]


def test_distal_subblock_snaps_to_containing_satellite(monkeypatch):
    starts, ends = [100, 1000], [500, 1400]
    names, monomers = [None, None], [171, 171]
    periodicities, hor_flags = [None, None], [False, False]
    link = DistalSatelliteLink(
        150, 220, 1050, 1200, "candidate_to_all", 0.8, 0.9, 0.9, 100
    )
    monkeypatch.setattr(
        cli,
        "detect_precise_boundaries",
        lambda **kwargs: (_ for _ in ()).throw(
            AssertionError("contained endpoint was revalidated")
        ),
    )

    [snapped] = cli.promote_unmatched_distal_candidates(
        [link],
        starts,
        ends,
        names,
        monomers,
        periodicities,
        hor_flags,
        fasta_file="input.fa",
        seq_id="chr1",
        seq_len=2000,
        window=100,
        k=21,
    )

    assert (snapped.start1, snapped.end1) == (100, 500)
    assert (snapped.start2, snapped.end2) == (1000, 1400)
    assert list(zip(starts, ends)) == [(100, 500), (1000, 1400)]


def test_verbose_boundary_merge_omits_raw_spans_and_appended_messages(capsys):
    result = cli.merge_shared_boundaries(
        intervals=[((0, 10), 3), ((100, 120), 4)],
        prefix=0,
        window=10,
        verbose=True,
    )

    assert result == [(0, 10, 3), (100, 120, 4)]
    output = capsys.readouterr().out
    assert "0 10 3" not in output
    assert "100 120 4" not in output
    assert "Appended new" not in output


def test_main_build_db_creates_output_files(monkeypatch, tmp_path, capsys):
    saved = []
    bed_df = pl.DataFrame(
        {"chrom": ["chr1"], "start": [1], "end": [5], "name": ["AlphaSat"]}
    )

    monkeypatch.setattr(sys, "argv", ["anianns", "build_db", "-f", "input.fa", "-b", "input.bed", "-d", str(tmp_path)])
    monkeypatch.setattr(cli, "read_bed_files", lambda files: [bed_df])
    monkeypatch.setattr(cli, "check_bed_vs_indexed_fasta", lambda bed, fasta: True)
    monkeypatch.setattr(
        cli,
        "extract_regions_by_name",
        lambda df, fasta, kmer, verbose: {"alphasat": {11, 22, 33}},
    )
    monkeypatch.setattr(
        cli,
        "save_kmer_sets_shared_k",
        lambda kmer_dict, k, output_path: saved.append((kmer_dict, k, output_path)),
    )

    cli.main()

    assert saved == [
        ({"alphasat": {11, 22, 33}}, 21, str(tmp_path / "alphasat.db"))
    ]
    assert "Building a k-mer database" in capsys.readouterr().out


@pytest.mark.parametrize(
    ("plot", "distal"),
    [(False, False), (True, False), (False, True), (True, True)],
)
def test_main_annotate_writes_bed_and_summary(
    monkeypatch, tmp_path, plot, distal
):
    monkeypatch.setattr(cli.os.path, "isfile", lambda path: True)
    argv = [
        "anianns",
        "annotate",
        "-f",
        "input.fa",
        "-d",
        str(tmp_path),
        "--band",
        "0.00004",
        "--quiet",
    ]
    if plot:
        argv.append("--plot")
    if distal:
        argv.append("--distal")
    monkeypatch.setattr(
        sys,
        "argv",
        argv,
    )
    monkeypatch.setattr(cli, "get_input_headers", lambda fasta: [("input.fa", ["chr1"])])
    monkeypatch.setattr(cli.pysam, "FastaFile", lambda path: FakeFastaFile({"chr1": "A" * 100}))
    monkeypatch.setattr(
        cli,
        "iter_hashed_fasta_bands",
        lambda fasta, seq_id, plan, **kwargs: [
            HashedBand(
                index=request.index,
                start=request.start,
                hashes=(4,) * request.hash_count,
            )
            for request in plan.requests()
        ],
    )
    monkeypatch.setattr(
        cli,
        "build_kmer_sets",
        lambda kmer_list, max_len, win, interval, sketch: (
            [__import__("numpy").array([1], dtype="int32")],
            [__import__("numpy").array([1], dtype="int32")],
        ),
    )
    monkeypatch.setattr(
        cli,
        "build_kmer_sets_multi",
        lambda kmer_list, window_configs, sketch: {
            window: (
                [__import__("numpy").array([1], dtype="int32")],
                [__import__("numpy").array([1], dtype="int32")],
            )
            for window in window_configs
        },
    )
    calls = {
        "dense": 0,
        "plot": 0,
        "highlights": [],
        "diagonal": [],
        "edges": [],
        "cmaps": [],
        "colorbars": [],
        "dpis": [],
        "white_below": [],
    }

    def fake_dense(*args):
        calls["dense"] += 1
        return __import__("numpy").array([[100.0]])

    monkeypatch.setattr(cli, "intersection_matrix_thresholded", fake_dense)
    monkeypatch.setattr(cli, "intersection_matrix", lambda *args: fake_dense())
    monkeypatch.setattr(
        cli,
        "intersection_matrix_with_threshold",
        lambda *args: (
            fake_dense(),
            __import__("numpy").array([[True]], dtype=bool),
        ),
    )
    monkeypatch.setattr(
        cli,
        "get_diagonal_span_from_sets",
        lambda ov, nov, win, k, identity, zero_tol: [((0, 20), 12)],
    )
    def fake_plot(matrix, **kwargs):
        calls["plot"] += 1
        calls["highlights"].append(kwargs.get("highlight_ranges"))
        calls["diagonal"].append(kwargs.get("diagonal_ranges"))
        calls["edges"].append(kwargs.get("edge_overlay"))
        calls["cmaps"].append(kwargs.get("cmap"))
        calls["colorbars"].append(kwargs.get("show_colorbar"))
        calls["dpis"].append(kwargs.get("dpi"))
        calls["white_below"].append(kwargs.get("white_below"))
        Path(kwargs["save_path"]).write_bytes(b"png")

    monkeypatch.setattr(cli, "plot_matrix", fake_plot)
    monkeypatch.setattr(cli, "detect_matrix_distal_links", lambda *args: [])
    monkeypatch.setattr(cli, "get_diagonal_span", lambda matrix, win, zero_tol: [((0, 20), 12)])
    monkeypatch.setattr(
        cli,
        "merge_shared_boundaries",
        lambda intervals, prefix, window, verbose=False: [(0, 20, 12)],
    )
    monkeypatch.setattr(
        cli,
        "report_borders",
        lambda **kwargs: ([2], [18], ["HSAT"], [171], [171], [False]),
    )

    cli.main()

    bed_path = tmp_path / "chr1.bed"
    csv_path = tmp_path / "chr1.csv"

    assert bed_path.exists()
    assert csv_path.exists()

    bed_text = bed_path.read_text()
    assert "chr1\t2\t18\tHSAT" in bed_text

    with csv_path.open() as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["name"] == "HSAT"
    assert rows[0]["intervals"] == "2-18"
    dsu_rows = list(
        csv.DictReader((tmp_path / "satellite_dsu.tsv").open(), delimiter="\t")
    )
    assert len(dsu_rows) == 1
    assert dsu_rows[0]["chrom"] == "chr1"
    assert dsu_rows[0]["start"] == "2"
    assert dsu_rows[0]["end"] == "18"
    assert dsu_rows[0]["component_size"] == "1"
    dsu_text = (tmp_path / "satellite_dsu.txt").read_text()
    assert "Total satellites: 1" in dsu_text
    assert "chr1:2-18 name=HSAT" in dsu_text
    assert calls["dense"] == (2 if plot or distal else 0)
    assert calls["plot"] == (4 if plot else 0)
    if plot and not distal:
        assert calls["highlights"] == [[], None, [], None]
    if plot:
        assert len(calls["edges"]) == 4
        assert all(calls["edges"][index] is not None for index in (0, 2))
        assert all(calls["edges"][index] is None for index in (1, 3))
        assert all(calls["diagonal"][index] for index in (0, 2))
        assert all(calls["diagonal"][index] is None for index in (1, 3))
        assert [calls["cmaps"][index] for index in (1, 3)] == [
            "spectral_11_r",
            "spectral_11_r",
        ]
        assert [calls["colorbars"][index] for index in (1, 3)] == [True, True]
        assert [calls["white_below"][index] for index in (1, 3)] == [86, 86]
        assert all(calls["dpis"][raw] > calls["dpis"][raw - 1] for raw in (1, 3))
    heatmaps = sorted((tmp_path / "matrix_plots").glob("*.png"))
    assert len(heatmaps) == (4 if plot else 0)
    assert len(list((tmp_path / "matrix_plots").glob("*_identity.png"))) == (
        2 if plot else 0
    )
    assert not list(tmp_path.glob("*_matrix.npy"))
    assert not list(tmp_path.glob("*_distal_*_matrix.npy"))
    assert (tmp_path / "chr1_distal_neighborhood.npz").exists() is distal
    links_path = tmp_path / "chr1_distal_links.bedpe"
    assert links_path.exists() is distal
    if links_path.exists():
        assert links_path.read_text().startswith("#chrom1\tstart1\tend1")


def test_distal_option_can_be_combined_with_plot():
    parser = cli.get_parser()
    args = parser.parse_args(
        ["annotate", "-f", "input.fa", "--distal", "--plot"]
    )
    assert args.distal is True
    assert args.plot is True


def test_binary_plot_downsampling_uses_max_pooling():
    matrix = __import__("numpy").zeros((6, 6), dtype=bool)
    matrix[1, 5] = True

    pooled, scale = cli.downsample_binary_matrix(matrix, max_pixels=3)

    assert scale == 2
    assert pooled.shape == (3, 3)
    assert pooled[0, 2]
    assert pooled.sum() == 1


def test_numeric_plot_downsampling_preserves_peak_identity():
    matrix = __import__("numpy").arange(36, dtype=float).reshape(6, 6)

    pooled, scale = cli.downsample_numeric_matrix(matrix, max_pixels=3)

    assert scale == 2
    assert pooled.shape == (3, 3)
    assert pooled[0, 0] == 7
    assert pooled[-1, -1] == 35


def test_seam_candidates_merge_fragments_across_band_boundary():
    candidates = [
        (80_000, 100_000, 12),
        (100_000, 135_000, 20),
        (300_000, 320_000, 10),
    ]

    cli.incorporate_seam_candidates(
        candidates,
        [(90_000, 120_000, 24)],
        window=1000,
    )

    assert candidates == [
        (80_000, 135_000, 32),
        (300_000, 320_000, 10),
    ]


def test_main_annotate_reports_missing_fasta(monkeypatch, tmp_path, capsys):
    missing_fasta = tmp_path / "sample_hap1_.fa"
    monkeypatch.setattr(
        sys,
        "argv",
        ["anianns", "annotate", "-f", str(missing_fasta)],
    )

    import pytest

    with pytest.raises(SystemExit) as exc_info:
        cli.main()

    captured = capsys.readouterr()
    assert exc_info.value.code == 1
    assert captured.err == f"[ERROR] FASTA file does not exist: {missing_fasta}\n"


def test_validate_ntrprism_range_start_must_be_less_than_end():
    # start == end
    assert cli.validate_ntrprism_range(100, 100, 1000) is not None
    # start > end
    assert cli.validate_ntrprism_range(500, 100, 1000) is not None
    # valid case: start < end within bounds
    assert cli.validate_ntrprism_range(100, 500, 1000) is None


def test_validate_ntrprism_range_out_of_bounds(monkeypatch, capsys):
    # Range where end exceeds sequence length
    assert cli.validate_ntrprism_range(0, 2000, 1000) is not None
    # Range where start is negative
    assert cli.validate_ntrprism_range(-1, 500, 1000) is not None

    # Confirm main() exits with error when range is out of bounds for a real fasta
    fasta_path = str(
        Path(__file__).resolve().parents[1] / "sample_sequences" / "sample_hap1.fa"
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "anianns",
            "ntrprism",
            "-f",
            fasta_path,
            "-s",
            "sample_hap1",
            "--range",
            "0",
            "99999999",  # beyond the 14 000 000 bp sequence
        ],
    )
    import pytest
    with pytest.raises(SystemExit) as exc_info:
        cli.main()
    assert exc_info.value.code == 1
    assert "[ERROR]" in capsys.readouterr().err
