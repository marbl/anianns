import mmh3
import numpy as np

from anianns.kmer_utils import (
    build_kmer_sets,
    build_kmer_sets_multi,
    calculate_hash_distances,
    convert_set_list_to_sorted_arrays,
    generate_kmers_from_fasta,
    generate_kmers_from_fasta_forward_only,
    generate_kmers_from_fasta_reverse_only,
    read_sequence_kmers_from_file,
    remove_ambiguous_bases,
)


def test_remove_ambiguous_bases_filters_known_homopolymers():
    keep = mmh3.hash("AAAA")
    drop = mmh3.hash("NNNN")

    cleaned = remove_ambiguous_bases([keep, drop], 4)

    assert keep in cleaned
    assert drop not in cleaned


def test_compiled_hash_distances_match_consecutive_occurrences():
    hashes = np.array([5, 7, 5, 5, 7, 9], dtype=np.int32)
    assert calculate_hash_distances(hashes).tolist() == [2, 1, 3]


def test_convert_set_list_to_sorted_arrays_preserves_sorted_values():
    arrays = convert_set_list_to_sorted_arrays([{3, 1}, {2}])
    assert np.array_equal(arrays[0], np.array([1, 3], dtype=np.int32))
    assert np.array_equal(arrays[1], np.array([2], dtype=np.int32))


def test_build_kmer_sets_creates_overlap_and_non_overlap_windows():
    overlap, non_overlap = build_kmer_sets(
        kmer_list=[0, 4, 8, 12, 16, 20],
        max_len=3,
        window=2,
        interval=1,
        prepend=[24],
    )

    assert [arr.tolist() for arr in overlap] == [[4, 8], [4, 8, 12, 16]]
    assert [arr.tolist() for arr in non_overlap] == [[4, 24], [8, 12]]


def test_build_kmer_sets_supports_denser_modulo_two_sketches():
    overlap, non_overlap = build_kmer_sets(
        kmer_list=[0, 2, 4, 6, 8, 10],
        max_len=3,
        window=2,
        interval=1,
        sketch=2,
    )

    assert [arr.tolist() for arr in overlap] == [[2, 4], [2, 4, 6, 8]]
    assert [arr.tolist() for arr in non_overlap] == [[2], [4, 6]]


def test_multi_window_sets_match_independent_prefix_builds():
    hashes = np.arange(0, 80, 4, dtype=np.int32)
    configs = {
        4: (14, 4, 2),
        8: (20, 3, 4),
    }

    combined = build_kmer_sets_multi(hashes, configs, sketch=4)

    for window, (hash_count, max_len, interval) in configs.items():
        independent = build_kmer_sets(
            hashes[:hash_count], max_len, window, interval, sketch=4
        )
        for combined_sets, independent_sets in zip(combined[window], independent):
            assert len(combined_sets) == len(independent_sets)
            assert all(
                np.array_equal(left, right)
                for left, right in zip(combined_sets, independent_sets)
            )


def test_generate_kmers_handles_reverse_complements_and_ambiguous_bases():
    seq = "ACTGN"
    kmers = list(generate_kmers_from_fasta(seq, k=4, quiet=True))

    first = mmh3.hash("ACTG", seed=42)
    rc = mmh3.hash("CAGT", seed=42)
    assert kmers == [min(first, rc), 0]


def test_forward_and_reverse_only_kmer_generators_hash_expected_strings():
    seq = "ACTG"
    assert list(generate_kmers_from_fasta_forward_only(seq, 4, True)) == [
        mmh3.hash("ACTG", seed=42)
    ]
    assert list(generate_kmers_from_fasta_reverse_only(seq, 4, True)) == [
        mmh3.hash("CAGT", seed=42)
    ]


def test_kmer_generators_emit_progress_when_not_quiet(capsys):
    seq = "A" * 80
    forward = list(generate_kmers_from_fasta_forward_only(seq, 4, False))
    reverse = list(generate_kmers_from_fasta_reverse_only(seq, 4, False))
    canonical = list(generate_kmers_from_fasta(seq, 4, False))

    assert len(forward) == len(reverse) == len(canonical) == 77
    output = capsys.readouterr().out
    assert "Hashing forward k-mers" in output
    assert "Hashing reverse k-mers" in output
    assert "Hashing k-mers" in output


def test_kmer_generators_handle_short_sequences_with_progress_enabled(capsys):
    seq = "ACTGA"

    forward = list(generate_kmers_from_fasta_forward_only(seq, 4, False))
    reverse = list(generate_kmers_from_fasta_reverse_only(seq, 4, False))
    canonical = list(generate_kmers_from_fasta(seq, 4, False))

    assert len(forward) == len(reverse) == len(canonical) == 2
    assert "Hashing k-mers" in capsys.readouterr().out


def test_kmer_generators_return_empty_for_sequences_shorter_than_k():
    assert list(generate_kmers_from_fasta_forward_only("ACT", 4, False)) == []
    assert list(generate_kmers_from_fasta_reverse_only("ACT", 4, False)) == []
    assert list(generate_kmers_from_fasta("ACT", 4, False)) == []


def test_read_sequence_kmers_from_file_uses_fasta_fetch(monkeypatch):
    class FakeFasta:
        def __init__(self, filename):
            self.filename = filename

        def fetch(self, seqid):
            assert seqid == "chr1"
            return "ACTGA"

    monkeypatch.setattr("anianns.kmer_utils.pysam.FastaFile", FakeFasta)
    monkeypatch.setattr("anianns.kmer_utils.generate_kmers_from_fasta", lambda seq, ksize, quiet: [1, 2, 3])

    assert read_sequence_kmers_from_file("fake.fa", "chr1", 4, True) == [[1, 2, 3]]
