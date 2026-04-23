from pathlib import Path

import pytest

from anianns.build_kmer_db import save_kmer_sets_shared_k
from anianns.classification import (
    calculate_distances,
    classify_kmers,
    load_all_kmer_dbs,
    top_n_frequent_distances,
)


def test_load_all_kmer_dbs_reads_multiple_databases(tmp_path):
    save_kmer_sets_shared_k({"alpha": {1, 2}}, 21, str(tmp_path / "alpha.db"))
    save_kmer_sets_shared_k({"beta": {3, 4}}, 19, str(tmp_path / "beta.db"))

    loaded = load_all_kmer_dbs(str(tmp_path))

    assert loaded["alpha"] == (21, {"alpha": {1, 2}})
    assert loaded["beta"] == (19, {"beta": {3, 4}})


def test_load_all_kmer_dbs_validates_path(tmp_path):
    missing = tmp_path / "missing"
    with pytest.raises(FileNotFoundError):
        load_all_kmer_dbs(str(missing))

    file_path = tmp_path / "file.txt"
    file_path.write_text("x")
    with pytest.raises(ValueError):
        load_all_kmer_dbs(str(file_path))

    with pytest.raises(ValueError):
        load_all_kmer_dbs(str(tmp_path))


def test_classify_kmers_returns_best_match_and_sorted_scores(capsys):
    best, results = classify_kmers(
        query_kmers={1, 2, 3, 4},
        kmer_supersets={
            "alpha": {1, 2, 3, 9},
            "beta": {1},
            "gamma": {7, 8},
        },
        min_overlap=2,
        verbose=True,
    )

    assert best == "alpha"
    assert [name for name, _, _ in results] == ["alpha", "beta", "gamma"]
    assert results[0][1] == 3
    assert "Top 3 classification matches" in capsys.readouterr().out


def test_classify_kmers_handles_empty_queries():
    best, results = classify_kmers(set(), {"alpha": {1, 2}})
    assert best == "Unclassified"
    assert results == []


def test_distance_helpers_report_repeated_spacing(capsys):
    distances = calculate_distances([5, 1, 5, 2, 5, 1])
    assert distances == [2, 2, 4]

    top = top_n_frequent_distances([2, 2, 4, 4, 4], n=2)
    assert top is None
    assert "Top 2 distances" in capsys.readouterr().out
