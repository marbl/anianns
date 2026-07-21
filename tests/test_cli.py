import csv
import sys
from pathlib import Path

import polars as pl
import pytest

from anianns import anianns as cli
from anianns.kmer_pipeline import HashedBand


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
    assert args.output_format == "bed"

    dense_args = parser.parse_args(
        ["annotate", "-f", "input.fa", "--sketch", "2"]
    )
    assert dense_args.sketch == 2
    assert args.band == 2.0


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


@pytest.mark.parametrize("plot", [False, True])
def test_main_annotate_writes_bed_and_summary(monkeypatch, tmp_path, plot):
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
        lambda fasta, seq_id, plan, cache_dir: [
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
        lambda kmer_list, max_len, win, interval, sketch: ([{1}], [{1}]),
    )
    monkeypatch.setattr(
        cli,
        "intersection_matrix_thresholded",
        lambda ov, nov, k, identity: __import__("numpy").array(
            [[True]], dtype="bool"
        ),
    )
    monkeypatch.setattr(
        cli,
        "intersection_matrix",
        lambda ov, nov, k: __import__("numpy").array([[100.0]]),
    )
    calls = {"combined": 0, "split": 0}

    def fake_combined(*args):
        calls["combined"] += 1
        return __import__("numpy").array([[100.0, 0.0], [0.0, 100.0]])

    def fake_split(matrix):
        calls["split"] += 1
        return matrix, matrix

    monkeypatch.setattr(
        cli,
        "intersection_matrix_inverted",
        fake_combined,
    )
    monkeypatch.setattr(cli, "split_diagonal_attached", fake_split)
    monkeypatch.setattr(cli, "plot_matrix", lambda matrix: None)
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
    assert calls == {
        "combined": 1 if plot else 0,
        "split": 1 if plot else 0,
    }


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
    assert "[ERROR]" in capsys.readouterr().out
