import csv
from datetime import datetime
import subprocess
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
    assert args.verbose is False
    assert args.quiet is False
    assert args.log is False
    assert args.cache_dir is None
    assert not hasattr(args, "no_cache")
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

    with pytest.raises(SystemExit):
        parser.parse_args(["annotate", "-f", "input.fa", "--extend"])

    logged_args = parser.parse_args(
        ["annotate", "-f", "input.fa", "--verbose", "--log"]
    )
    assert logged_args.verbose is True
    assert logged_args.log is True

    cache_args = parser.parse_args(
        ["annotate", "-f", "input.fa", "--cache-dir", "cache"]
    )
    assert cache_args.cache_dir == "cache"
    with pytest.raises(SystemExit):
        parser.parse_args(["annotate", "-f", "input.fa", "--no-cache"])

    multi_args = parser.parse_args(
        ["annotate", "-f", "input.fa", "-w", "5000"]
    )
    assert cli.derive_window_sizes(multi_args.window) == (2500, 5000, 10_000)

    with pytest.raises(SystemExit):
        parser.parse_args(
            ["annotate", "-f", "input.fa", "-w", "1000", "2000"]
        )


def test_format_sequence_size_uses_scaled_decimal_units():
    assert cli.format_sequence_size(500) == "0.5kb"
    assert cli.format_sequence_size(750_000) == "750kb"
    assert cli.format_sequence_size(12_345_678) == "12.3mb"
    assert cli.format_sequence_size(100_000_000) == "100mb"
    assert cli.format_sequence_size(1_000_000_000) == "1gb"
    assert cli.format_sequence_size(1_250_000_000) == "1.25gb"
    with pytest.raises(ValueError, match="cannot be negative"):
        cli.format_sequence_size(-1)


def test_annotate_log_filename_contains_run_date_and_time():
    run_time = datetime(2026, 8, 4, 14, 37, 52)

    assert cli.annotate_log_filename(run_time) == (
        "anianns_annotation_log_2026-08-04_14-37-52.txt"
    )


def test_cache_hit_is_announced_immediately_before_matrix_creation(
    tmp_path, capsys
):
    cache_dir = tmp_path / "hashes"

    cli.announce_matrix_creation("chr2", 100_000_000, cache_dir)

    assert capsys.readouterr().out == (
        f"Found hashes in {cache_dir}\n"
        "Creating an ANI matrix for chr2 (100mb):\n\n"
    )


def test_matrix_creation_without_cache_has_no_hash_notice(capsys):
    cli.announce_matrix_creation("chr2", 100_000_000)

    assert capsys.readouterr().out == (
        "Creating an ANI matrix for chr2 (100mb):\n\n"
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


def test_rejects_unmatched_distal_candidate_when_ntr_fails_and_link_is_weak(
    monkeypatch,
):
    starts, ends = [100], [500]
    names, monomers, periodicities, hor_flags = [None], [171], [None], [False]
    link = DistalSatelliteLink(
        100, 500, 1000, 1400, "candidate_to_all", 0.4, 0.9, 0.9, 100
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


def test_strong_distal_link_cannot_override_ntr_rejection(monkeypatch):
    starts, ends = [100], [500]
    names, monomers, periodicities, hor_flags = ["telomere"], [6], [None], [False]
    link = DistalSatelliteLink(
        100, 500, 1000, 1400, "component", 0.55, 1.0, 1.0, 100
    )
    observed = {}

    def fake_boundaries(**kwargs):
        observed.update(kwargs)
        return None

    monkeypatch.setattr(cli, "detect_precise_boundaries", fake_boundaries)

    promoted = cli.promote_unmatched_distal_candidates(
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

    assert observed["coordinates"] == (1000, 1400)
    assert "allow_ntr_rejection" not in observed
    assert "fallback_prism_result" not in observed
    assert list(zip(starts, ends)) == [(100, 500)]
    assert monomers == [6]
    assert promoted == []


def test_ntr_rejected_endpoint_is_not_retried_for_stronger_link(monkeypatch):
    starts, ends = [100], [500]
    names, monomers, periodicities, hor_flags = [None], [6], [None], [False]
    weak = DistalSatelliteLink(
        100, 500, 1000, 1400, "candidate_to_all", 0.4, 0.9, 0.9, 100
    )
    strong = DistalSatelliteLink(
        100, 500, 1000, 1400, "component", 0.55, 1.0, 1.0, 100
    )
    attempts = []

    def fake_boundaries(**kwargs):
        attempts.append(kwargs["coordinates"])
        return None

    monkeypatch.setattr(cli, "detect_precise_boundaries", fake_boundaries)

    promoted = cli.promote_unmatched_distal_candidates(
        [weak, strong],
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

    assert attempts == [(1000, 1400)]
    assert list(zip(starts, ends)) == [(100, 500)]
    assert promoted == []


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


def test_resolve_overlapping_satellites_splits_unrelated_boundary_estimates():
    starts = [687654, 666461]
    ends = [697156, 688000]
    names = [None, None]
    monomers = [5, 37]
    periodicities = [None, None]
    hor_flags = [False, False]

    cli.resolve_overlapping_satellites(
        starts,
        ends,
        names,
        monomers,
        periodicities,
        hor_flags,
    )

    assert list(zip(starts, ends)) == [(666461, 687827), (687827, 697156)]
    assert monomers == [37, 5]
    assert all(left_end <= right_start for left_end, right_start in zip(ends, starts[1:]))


def test_resolve_overlapping_satellites_discards_contained_call():
    starts = [100, 200, 1000]
    ends = [900, 400, 1200]
    names = ["outer", "contained", "next"]
    monomers = [37, 5, 171]
    periodicities = [None, None, None]
    hor_flags = [False, False, False]

    cli.resolve_overlapping_satellites(
        starts,
        ends,
        names,
        monomers,
        periodicities,
        hor_flags,
    )

    assert list(zip(starts, ends)) == [(100, 900), (1000, 1200)]
    assert names == ["outer", "next"]
    assert monomers == [37, 171]


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
    ("plot", "distal", "output_format", "verbose"),
    [
        (False, False, "bed", True),
        (True, False, "bed", False),
        (False, True, "bed", False),
        (True, True, "bed", False),
        (False, False, "gtf", False),
    ],
)
def test_main_annotate_writes_requested_format_and_summary(
    monkeypatch,
    tmp_path,
    capsys,
    plot,
    distal,
    output_format,
    verbose,
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
        "--output-format",
        output_format,
    ]
    argv.append("--log" if verbose else "--quiet")
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
        "legends_outside": [],
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
        calls["legends_outside"].append(kwargs.get("legend_outside", False))
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
    expected_plot_calls = 4 if plot else 0

    annotation_path = tmp_path / f"chr1.{output_format}"
    csv_path = tmp_path / "chr1.csv"

    assert annotation_path.exists()
    assert csv_path.exists()

    annotation_text = annotation_path.read_text()
    if output_format == "bed":
        assert "chr1\t2\t18\tHSAT" in annotation_text
        annotation_fields = annotation_text.splitlines()[1].split("\t")
        assert annotation_fields[-1] != "0,0,0"
        assert len(annotation_fields[-1].split(",")) == 3
    else:
        assert annotation_text.startswith(
            "chr1\tAniAnns\ttandem_repeat\t3\t18\t.\t.\t.\t"
        )
        assert 'repeat_name "HSAT";' in annotation_text

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
    assert calls["plot"] == expected_plot_calls
    if plot and not distal:
        assert calls["highlights"] == [[], None, [], None]
    if plot:
        raw_indexes = tuple(range(0, expected_plot_calls, 2))
        identity_indexes = tuple(range(1, expected_plot_calls, 2))
        assert len(calls["edges"]) == expected_plot_calls
        assert all(calls["edges"][index] is not None for index in raw_indexes)
        assert all(calls["edges"][index] is None for index in identity_indexes)
        assert all(calls["diagonal"][index] for index in raw_indexes)
        assert all(calls["diagonal"][index] is None for index in identity_indexes)
        assert all(
            calls["cmaps"][index] == "spectral_11_r"
            for index in identity_indexes
        )
        assert all(calls["colorbars"][index] for index in identity_indexes)
        assert all(calls["white_below"][index] == 86 for index in identity_indexes)
        assert calls["legends_outside"] == [True, False] * len(raw_indexes)
        assert all(
            calls["dpis"][identity] > calls["dpis"][identity - 1]
            for identity in identity_indexes
        )
    heatmaps = sorted((tmp_path / "matrix_plots").glob("*.png"))
    pair_heatmaps = sorted((tmp_path / "matrix_pairs").glob("*.png"))
    expected_band_plots = 4 if plot else 0
    assert len(heatmaps) == expected_band_plots
    assert pair_heatmaps == []
    assert len(list((tmp_path / "matrix_plots").glob("*_identity.png"))) == (
        expected_band_plots // 2
    )
    assert not list(tmp_path.glob("*_matrix.npy"))
    assert not list(tmp_path.glob("*_distal_*_matrix.npy"))
    assert (
        tmp_path / "chr1_distal_neighborhood.npz"
    ).exists() is distal
    links_path = tmp_path / "chr1_distal_links.bedpe"
    assert links_path.exists() is distal
    if links_path.exists():
        assert links_path.read_text().startswith("#chrom1\tstart1\tend1")
    captured = capsys.readouterr()
    if verbose:
        assert " potential candidates" in captured.out
        assert "Matrix 1/" not in captured.out
        log_paths = list(tmp_path.glob("anianns_annotation_log_*.txt"))
        assert len(log_paths) == 1
        log_path = log_paths[0]
        assert log_path.exists()
        log_text = log_path.read_text()
        assert " potential candidates" in log_text
        assert "Matrix 1/" not in log_text
        assert cli.VERSION in log_text
    else:
        assert captured.out == ""
        assert captured.err == ""


def test_quiet_and_log_conflict_is_silent(monkeypatch, capsys):
    monkeypatch.setattr(
        sys,
        "argv",
        ["anianns", "annotate", "-f", "input.fa", "--quiet", "--log"],
    )

    assert cli.main() == 2
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


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


def test_quiet_mode_is_silent_for_entire_cli_process(tmp_path):
    missing_fasta = tmp_path / "missing.fa"

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "anianns",
            "annotate",
            "-f",
            str(missing_fasta),
            "--quiet",
        ],
        capture_output=True,
        check=False,
    )

    assert result.returncode == 1
    assert result.stdout == b""
    assert result.stderr == b""


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
