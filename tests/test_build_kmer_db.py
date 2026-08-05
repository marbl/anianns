import gzip
import struct

from anianns.build_kmer_db import load_kmer_sets_shared_k, save_kmer_sets_shared_k


def test_save_and_load_kmer_sets_round_trip(tmp_path):
    path = tmp_path / "satellites.db"
    original = {"alpha": {1, -2, 3}, "beta": {42}}

    save_kmer_sets_shared_k(original, k=21, output_path=str(path))
    k, loaded = load_kmer_sets_shared_k(str(path))

    assert k == 21
    assert loaded == original


def test_load_kmer_sets_shared_k_returns_empty_result_on_truncated_file(
    tmp_path, capsys
):
    path = tmp_path / "broken.db"
    with gzip.open(path, "wb") as handle:
        handle.write(struct.pack("<B", 21))
        handle.write(struct.pack("<H", 1))

    k, loaded = load_kmer_sets_shared_k(str(path))

    assert k is None
    assert loaded == {}
    assert "Failed to load file" in capsys.readouterr().out
