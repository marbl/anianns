import sys

import numpy as np
import pysam
import pytest

from anianns import anianns as cli
from anianns.ntrprism import (
    SpacingPeak,
    analyze_kmer_spacings,
    format_ascii_histogram,
    format_spacing_table,
    merge_spacing_counts,
)


def write_indexed_fasta(tmp_path, sequence="ACGT" * 100):
    fasta_path = tmp_path / "input.fa"
    fasta_path.write_text(f">chrTest\n{sequence}\n")
    pysam.faidx(str(fasta_path))
    return fasta_path


def test_merge_spacing_counts_stacks_neighboring_values():
    peaks = merge_spacing_counts({170: 5, 171: 7, 173: 4, 200: 2}, tolerance=1)

    assert peaks[0].spacing == 171
    assert peaks[0].minimum == 170
    assert peaks[0].maximum == 171
    assert peaks[0].count == 12
    assert [(peak.spacing, peak.count) for peak in peaks[1:]] == [
        (173, 4),
        (200, 2),
    ]


def test_analyze_kmer_spacings_finds_tandem_period():
    peaks, total = analyze_kmer_spacings("ACGT" * 20, kmer=3, merge_distance=1)

    assert total == 74
    assert peaks[0].spacing == 4
    assert peaks[0].count == 74


def test_analyze_kmer_spacings_defaults_to_k21(monkeypatch):
    observed = {}

    def capture_hashes(sequence, kmer):
        observed["kmer"] = kmer
        return np.arange(len(sequence) - kmer + 1, dtype=np.int32)

    monkeypatch.setattr("anianns.ntrprism.forward_kmer_hashes", capture_hashes)

    analyze_kmer_spacings("ACGT" * 20)

    assert observed["kmer"] == 21


def test_ntrprism_parser_defaults():
    args = cli.get_parser().parse_args(
        ["ntrprism", "-f", "input.fa", "-s", "chr1", "-r", "10", "100"]
    )

    assert args.fasta == "input.fa"
    assert args.seq_id == "chr1"
    assert args.range == [10, 100]
    assert args.kmer == 21
    assert args.merge_distance == 1
    assert args.save is False


def test_annotate_parser_defaults_to_k21():
    args = cli.get_parser().parse_args(["annotate", "-f", "input.fa"])

    assert args.kmer == 21


def test_ntrprism_parser_accepts_custom_kmer_and_save():
    args = cli.get_parser().parse_args(
        [
            "ntrprism",
            "-f",
            "input.fa",
            "-s",
            "chr1",
            "-r",
            "10",
            "100",
            "-k",
            "11",
            "--save",
        ]
    )

    assert args.kmer == 11
    assert args.save is True


def test_terminal_spacing_table_and_ascii_histogram_use_interval_percentage():
    peaks = [
        SpacingPeak(spacing=2241, minimum=2240, maximum=2242, count=25),
        SpacingPeak(spacing=2, minimum=1, maximum=3, count=5),
    ]

    table = format_spacing_table(peaks, interval_length=100, top_n=10)
    histogram = format_ascii_histogram(peaks, interval_length=100, top_n=10)

    assert "2,241" in table
    assert "25.000%" in table
    assert "2,241 bp | ################################################" in histogram
    assert "2 bp | ##########" in histogram
    assert "5.000%" in histogram


def test_main_ntrprism_reports_missing_fasta(monkeypatch, tmp_path, capsys):
    missing = tmp_path / "missing.fa"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "anianns",
            "ntrprism",
            "-f",
            str(missing),
            "-s",
            "chr1",
            "-r",
            "0",
            "100",
        ],
    )

    with pytest.raises(SystemExit) as exc_info:
        cli.main()

    assert exc_info.value.code == 1
    assert (
        f"[ERROR] FASTA file does not exist: {missing}"
        in capsys.readouterr().err
    )


def test_main_ntrprism_reports_missing_sequence(monkeypatch, tmp_path, capsys):
    fasta_path = write_indexed_fasta(tmp_path)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "anianns",
            "ntrprism",
            "-f",
            str(fasta_path),
            "-s",
            "missing",
            "-r",
            "0",
            "100",
        ],
    )

    with pytest.raises(SystemExit) as exc_info:
        cli.main()

    assert exc_info.value.code == 1
    assert "Sequence ID 'missing' was not found" in capsys.readouterr().err


def test_main_ntrprism_prints_report_without_saving(monkeypatch, tmp_path, capsys):
    fasta_path = write_indexed_fasta(tmp_path)
    output_dir = tmp_path / "output"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "anianns",
            "ntrprism",
            "-f",
            str(fasta_path),
            "-s",
            "chrTest",
            "-r",
            "0",
            "400",
            "-k",
            "3",
            "-d",
            str(output_dir),
        ],
    )

    result = cli.main()
    output = capsys.readouterr().out

    assert result == 0
    assert "(400 bp, k=3)" in output
    assert "Top k-mer spacing distances" in output
    assert "98.500%" in output
    assert "ASCII histogram" in output
    assert "4 bp |" in output
    assert not output_dir.exists()


def test_main_ntrprism_writes_files_only_with_save(monkeypatch, tmp_path):
    fasta_path = write_indexed_fasta(tmp_path)
    output_dir = tmp_path / "output"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "anianns",
            "ntrprism",
            "-f",
            str(fasta_path),
            "-s",
            "chrTest",
            "-r",
            "0",
            "400",
            "-k",
            "3",
            "-d",
            str(output_dir),
            "--save",
            "--quiet",
        ],
    )

    result = cli.main()

    report_text = (output_dir / "chrTest_0_400_ntrprism.txt").read_text()
    assert result == 0
    assert "Interval_percentage\tRepeated_spacing_fraction" in report_text
    assert "Kmer_length\t3" in report_text
    assert "1\t4\t4\t394\t98.500000" in report_text
    assert "ASCII histogram" in report_text
    assert "4 bp |" in report_text
    assert (output_dir / "chrTest_0_400_ntrprism.png").stat().st_size > 0
