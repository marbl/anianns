import csv
import json

import polars as pl
import pytest

from anianns.general_utils import (
    add_prefix_to_tuples,
    calculate_distances,
    check_bed_vs_indexed_fasta,
    clean_genomic_ticks,
    convert_dataframe_format,
    define_bounds,
    extract_histograms_by_name,
    extract_region,
    extract_regions_by_name,
    get_fasta_indexed_chroms,
    get_input_headers,
    merge_close_values,
    plot_matrix,
    read_bed_files,
    validate_json,
    write_summary_file,
)


def test_clean_genomic_ticks_uses_band_coordinates_and_readable_steps():
    assert clean_genomic_ticks(2.0, 4.0).tolist() == [2.0, 2.5, 3.0, 3.5, 4.0]


@pytest.fixture
def bed_like_df():
    return pl.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "start": [10, 20],
            "end": [15, 25],
            "name": ["alpha", "beta"],
            "score": [0, 1],
            "strand": [".", "+"],
            "thickStart": [10, 20],
            "thickEnd": [15, 25],
            "itemRgb": ["255,0,0", "0,0,0"],
        }
    )


def test_add_prefix_to_tuples_and_calculate_distances():
    assert add_prefix_to_tuples([(1, 5), (10, 15)], band_height=100, w=3) == [
        (201, 205),
        (210, 215),
    ]
    assert calculate_distances([1, 3, 1, 1, 3]) == [2, 1, 3]


def test_check_bed_vs_indexed_fasta(monkeypatch, capsys):
    monkeypatch.setattr(
        "anianns.general_utils.get_fasta_indexed_chroms",
        lambda paths: {"chr1", "chr2"},
    )
    dfs = [pl.DataFrame({"chrom": ["chr1"], "start": [1], "end": [2]})]
    assert check_bed_vs_indexed_fasta(dfs, ["fake.fa"]) is True

    missing = [pl.DataFrame({"chrom": ["chr3"], "start": [1], "end": [2]})]
    assert check_bed_vs_indexed_fasta(missing, ["fake.fa"]) is False
    assert "Input fastas missing" in capsys.readouterr().out


def test_convert_dataframe_format_supports_multiple_targets(bed_like_df):
    gtf = convert_dataframe_format(bed_like_df, "gtf")
    first_gtf = gtf.splitlines()[0].split("\t")
    assert first_gtf[:5] == ["chr1", "AniAnns", "tandem_repeat", "11", "15"]
    assert first_gtf[5:8] == [".", ".", "."]
    assert 'gene_id "anianns_000001";' in first_gtf[8]
    assert 'repeat_name "alpha";' in first_gtf[8]
    assert 'monomer_length "0";' in first_gtf[8]

    gff = convert_dataframe_format(bed_like_df, "gff")
    assert gff.startswith("##gff-version 3\n")
    assert "ID=anianns_000002;Name=beta;monomer_length=1" in gff

    csv_output = convert_dataframe_format(bed_like_df, "csv")
    assert "chrom,start,end,name" in csv_output

    tsv_output = convert_dataframe_format(bed_like_df, "tsv")
    assert "chrom\tstart\tend\tname" in tsv_output

    json_output = convert_dataframe_format(bed_like_df, "json")
    assert '"chrom":"chr1"' in json_output

    with pytest.raises(ValueError):
        convert_dataframe_format(bed_like_df, "xlsx")


def test_convert_dataframe_format_accepts_final_bed_chrom_column():
    output = convert_dataframe_format(
        pl.DataFrame(
            {
                "#chrom": ["chr2_hap1"],
                "start": [666461],
                "end": [687827],
                "name": [None],
                "score": [37],
                "strand": ["."],
            }
        ),
        "gtf",
    )

    fields = output.rstrip("\n").split("\t")
    assert fields[:5] == [
        "chr2_hap1",
        "AniAnns",
        "tandem_repeat",
        "666462",
        "687827",
    ]
    assert 'repeat_name "Unclassified Repeat";' in fields[8]


def test_define_bounds_parses_supported_identifiers():
    assert define_bounds("chrY:50-3000") == ("chrY", 50, 3000)
    assert define_bounds("HG002_chr13_MATERNAL:1-4000000:1000000-3000000") == (
        "HG002_chr13_MATERNAL",
        1000000,
        3000000,
    )
    assert define_bounds("chr1") is None


def test_merge_close_values_groups_neighbors_by_tolerance():
    merged = merge_close_values([(100, 2), (102, 5), (200, 1), (201, 4)], tolerance=3)
    assert merged == [(102, 7), (201, 5)]


def test_extract_region_returns_sequence_and_handles_errors(monkeypatch, capsys):
    class FakeFasta:
        def fetch(self, chrom, start, end):
            assert start == 1
            return f"{chrom}:{start}-{end}"

        def close(self):
            return None

    monkeypatch.setattr(
        "anianns.general_utils.pysam.FastaFile", lambda path: FakeFasta()
    )
    assert extract_region("fake.fa", "chr1", 0, 10) == "chr1:1-10"

    def raise_error(path):
        raise RuntimeError("bad fasta")

    monkeypatch.setattr("anianns.general_utils.pysam.FastaFile", raise_error)
    assert extract_region("fake.fa", "chr1", 1, 10) is None
    assert "Error fetching region" in capsys.readouterr().out


def test_extract_region_reuses_an_open_fasta_without_closing_it():
    class OpenFasta:
        def __init__(self):
            self.closed = False

        def fetch(self, chrom, start, end):
            return f"{chrom}:{start}-{end}"

        def close(self):
            self.closed = True

    fasta = OpenFasta()
    assert extract_region(fasta, "chr1", 1, 10) == "chr1:1-10"
    assert fasta.closed is False


def test_extract_regions_and_histograms_by_name(monkeypatch):
    df = pl.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "start": [1, 5],
            "end": [6, 10],
            "name": ["Alpha", "Alpha"],
        }
    )
    calls = []

    def fake_extract_region(fasta, chrom, start, end):
        calls.append((fasta, chrom, start, end))
        return "ACTGA"

    monkeypatch.setattr("anianns.general_utils.extract_region", fake_extract_region)
    monkeypatch.setattr(
        "anianns.general_utils.generate_kmers_from_fasta",
        lambda seq, k, quiet: [1, 2, 3],
    )
    monkeypatch.setattr(
        "anianns.general_utils.generate_kmers_from_fasta_forward_only",
        lambda seq, k, quiet: iter([9, 1, 9, 1, 9]),
    )
    monkeypatch.setattr(
        "anianns.general_utils.top_n_frequent_distances", lambda values, n: [(2, 3)]
    )

    regions = extract_regions_by_name(df, ["one.fa", "two.fa"], k=4, verbose=False)
    histograms = extract_histograms_by_name(df, "one.fa", k=4, verbose=False)

    assert regions == {"alpha": {1, 2, 3}}
    assert histograms == {"alpha": [[(2, 3)], [(2, 3)]]}
    assert calls[0] == ("one.fa", "chr1", 1, 6)


def test_get_fasta_indexed_chroms_and_headers(monkeypatch, capsys):
    class FakeFasta:
        def __init__(self, path):
            if path == "bad.fa":
                raise OSError("missing")
            self.references = ["chr1", "chr2"]

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    monkeypatch.setattr("anianns.general_utils.pysam.FastaFile", FakeFasta)

    assert get_fasta_indexed_chroms(["good.fa", "bad.fa"]) == {"chr1", "chr2"}
    assert get_input_headers(["good.fa", "bad.fa"]) == [("good.fa", ["chr1", "chr2"])]
    assert "Failed to open or index FASTA file" in capsys.readouterr().out


def test_plot_matrix_validates_inputs(monkeypatch, tmp_path):
    monkeypatch.setattr("anianns.general_utils.plt.show", lambda: None)

    with pytest.raises(TypeError):
        plot_matrix([[1, 2]])

    with pytest.raises(ValueError):
        plot_matrix(pl.Series("x", [1, 2]).to_numpy())

    with pytest.raises(ValueError, match="overlay must match"):
        plot_matrix(
            pl.DataFrame([[1, 2], [3, 4]]).to_numpy(),
            edge_overlay=pl.DataFrame([[True]]).to_numpy(),
        )

    with pytest.raises(ValueError, match="requires a Sobel"):
        plot_matrix(
            pl.DataFrame([[1, 2], [3, 4]]).to_numpy(),
            edge_only=True,
        )

    plot_matrix(pl.DataFrame([[1, 2], [3, 4]]).to_numpy(), show_colorbar=False)

    save_path = tmp_path / "heatmap.png"
    monkeypatch.setattr(
        "anianns.general_utils.plt.show",
        lambda: (_ for _ in ()).throw(AssertionError("interactive display used")),
    )
    plot_matrix(
        pl.DataFrame([[86, 90], [95, 100]]).to_numpy(),
        edge_overlay=pl.DataFrame([[True, False], [False, True]]).to_numpy(),
        edge_only=True,
        diagonal_ranges=[(0, 1)],
        highlight_ranges=[(1, 2, 0, 1)],
        show_colorbar=False,
        dpi=36,
        figsize=(2, 2),
        save_path=save_path,
    )
    assert save_path.exists()


def test_plot_matrix_can_place_overlay_legend_outside_axes(monkeypatch, tmp_path):
    captured = {}
    monkeypatch.setattr(
        "anianns.general_utils.plt.close",
        lambda figure: captured.setdefault("figure", figure),
    )

    plot_matrix(
        pl.DataFrame([[0, 1], [1, 0]]).to_numpy(),
        edge_overlay=pl.DataFrame([[True, False], [False, True]]).to_numpy(),
        edge_only=True,
        diagonal_ranges=[(0, 1)],
        highlight_ranges=[(1, 2, 0, 1)],
        legend_outside=True,
        show_colorbar=False,
        save_path=tmp_path / "outside-legend.png",
    )

    figure = captured["figure"]
    assert figure.axes[0].get_legend() is None
    assert len(figure.legends) == 1
    labels = [text.get_text() for text in figure.legends[0].get_texts()]
    assert labels == ["Sobel edges", "Detected satellite", "Distal link"]


def test_read_bed_files_skips_headers_and_casts_coordinates(tmp_path):
    bed_path = tmp_path / "example.bed"
    bed_path.write_text(
        "track name=test\n" "# ignored\n" "chr1\t10\t20\talpha\n" "chr2\t30\t40\tbeta\n"
    )

    [df] = read_bed_files(str(bed_path))

    assert df.columns == ["chrom", "start", "end", "name"]
    assert df["start"].dtype == pl.Int64
    assert df["name"].to_list() == ["alpha", "beta"]


def test_validate_json_checks_existence_syntax_and_required_keys(tmp_path, capsys):
    missing = tmp_path / "missing.json"
    assert validate_json(str(missing)) is False

    invalid = tmp_path / "invalid.json"
    invalid.write_text("{bad json")
    assert validate_json(str(invalid)) is False

    valid = tmp_path / "valid.json"
    valid.write_text(json.dumps({"alpha": 1, "beta": 2}))
    assert validate_json(str(valid), required_keys=["alpha"]) is True
    assert validate_json(str(valid), required_keys=["gamma"]) is False

    output = capsys.readouterr().out
    assert "File does not exist" in output
    assert "Invalid JSON syntax" in output
    assert "Missing required keys" in output


def test_write_summary_file_groups_named_rows_and_splits_unknowns(tmp_path):
    out_csv = tmp_path / "summary.csv"
    write_summary_file(
        (
            [10, 30, 50],
            [20, 40, 60],
            ["alpha", "alpha", "Unknown"],
            [171, 171, 180],
            [171, 171, 180],
            [False, False, True],
        ),
        str(out_csv),
    )

    with out_csv.open() as handle:
        rows = list(csv.DictReader(handle))

    assert rows[0]["name"] == "alpha"
    assert rows[0]["intervals"] == "10-20;30-40"
    assert rows[1]["name"] == "Unclassified Repeat"


def test_write_summary_file_requires_aligned_lengths(tmp_path):
    with pytest.raises(ValueError):
        write_summary_file(([1], [2], ["x"], [3], [4], []), str(tmp_path / "x.csv"))
