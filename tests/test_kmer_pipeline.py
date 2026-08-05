import numpy as np
import pysam

from anianns import kmer_pipeline
from anianns.kmer_pipeline import (
    HashedBand,
    SequenceBandPlan,
    canonical_kmer_hashes,
    forward_kmer_hashes,
    iter_hashed_fasta_bands,
    load_cached_sequence_hashes,
)
from anianns.kmer_utils import generate_kmers_from_fasta
from anianns.kmer_utils import generate_kmers_from_fasta_forward_only


def test_canonical_kmer_hashes_match_existing_generator_exactly():
    sequences = [
        "ACTGactgACTG",
        "AAARYMKSWNtttACTG",
        "GATTACAGATTACA",
        "ACT",
    ]

    for sequence in sequences:
        expected = np.array(
            list(generate_kmers_from_fasta(sequence, 4, True)), dtype=np.int32
        )
        assert np.array_equal(canonical_kmer_hashes(sequence, 4), expected)


def test_canonical_kmer_hashes_match_reference_on_randomized_sequences():
    rng = np.random.default_rng(2026)
    alphabet = np.array(list("ACGTrymkswhbvdn-"))
    for length in (1, 4, 21, 64, 257):
        sequence = "".join(rng.choice(alphabet, size=length).tolist())
        for kmer in (1, 4, 21):
            expected = np.array(
                list(generate_kmers_from_fasta(sequence, kmer, True)),
                dtype=np.int32,
            )
            assert np.array_equal(canonical_kmer_hashes(sequence, kmer), expected)


def test_batch_forward_hashes_match_existing_generator_exactly():
    for sequence in ("ACTGactgACTG", "AAARYMKSWNtttACTG", "ACT"):
        expected = np.array(
            list(generate_kmers_from_fasta_forward_only(sequence, 6, True)),
            dtype=np.int32,
        )
        assert np.array_equal(forward_kmer_hashes(sequence, 6), expected)


def test_sequence_band_plan_covers_tail_and_shares_multi_window_hashes():
    plan = SequenceBandPlan(
        sequence_length=105,
        kmer=5,
        band_height=40,
        windows=(10, 20, 10),
    )

    assert plan.total_kmers == 101
    assert plan.band_count == 3
    assert plan.windows == (10, 20)
    assert plan.max_interval == 10
    assert [(r.index, r.start, r.hash_count) for r in plan.requests()] == [
        (0, 0, 50),
        (1, 40, 50),
        (2, 80, 21),
    ]

    shared_hashes = np.arange(50, dtype=np.int32)
    shared = HashedBand(index=0, start=0, hashes=shared_hashes)
    window_10_hashes = plan.hashes_for_window(shared, 10)
    window_20_hashes = plan.hashes_for_window(shared, 20)
    assert len(window_10_hashes) == 45
    assert len(window_20_hashes) == 50
    assert np.shares_memory(window_10_hashes, shared_hashes)
    assert np.shares_memory(window_20_hashes, shared_hashes)

    tail = HashedBand(index=2, start=80, hashes=tuple(range(21)))
    assert len(plan.hashes_for_window(tail, 10)) == 21
    assert len(plan.hashes_for_window(tail, 20)) == 21


def _indexed_fasta(tmp_path):
    fasta_path = tmp_path / "input.fa"
    sequence = "ACTGNNacgtRYMKSWACTGACTGACGT"
    fasta_path.write_text(f">chr1\n{sequence}\n")
    pysam.faidx(str(fasta_path))
    return fasta_path, sequence


def _assert_pipeline_matches_full_sequence(fasta_path, sequence, use_process):
    plan = SequenceBandPlan(
        sequence_length=len(sequence),
        kmer=4,
        band_height=10,
        windows=(6,),
    )
    expected = np.array(
        list(generate_kmers_from_fasta(sequence, 4, True)), dtype=np.int32
    )
    bands = list(
        iter_hashed_fasta_bands(str(fasta_path), "chr1", plan, use_process=use_process)
    )

    assert [band.index for band in bands] == list(range(plan.band_count))
    for request, band in zip(plan.requests(), bands):
        assert band.start == request.start
        assert np.array_equal(
            band.hashes,
            expected[request.start : request.start + request.hash_count],
        )


def test_synchronous_band_pipeline_matches_full_sequence_hashes(tmp_path):
    fasta_path, sequence = _indexed_fasta(tmp_path)
    _assert_pipeline_matches_full_sequence(fasta_path, sequence, use_process=False)


def test_hash_worker_applies_its_thread_limit(tmp_path, monkeypatch):
    fasta_path, sequence = _indexed_fasta(tmp_path)
    plan = SequenceBandPlan(len(sequence), 4, 10, (6,))
    request = next(plan.requests())
    configured = []
    monkeypatch.setattr(kmer_pipeline, "set_num_threads", configured.append)

    band = kmer_pipeline._hash_fasta_band(
        str(fasta_path), "chr1", len(sequence), 4, request, 1
    )

    assert configured == [1]
    assert band.index == request.index
    assert len(band.hashes) == request.hash_count


def test_process_band_pipeline_matches_full_sequence_hashes(tmp_path):
    fasta_path, sequence = _indexed_fasta(tmp_path)
    _assert_pipeline_matches_full_sequence(fasta_path, sequence, use_process=True)


def test_disk_hash_cache_is_reused_across_band_and_window_plans(tmp_path, monkeypatch):
    fasta_path, sequence = _indexed_fasta(tmp_path)
    cache_dir = tmp_path / "cache"
    expected = canonical_kmer_hashes(sequence, 4)
    calls = 0
    original_hash_band = kmer_pipeline._hash_fasta_band

    def counted_hash_band(*args):
        nonlocal calls
        calls += 1
        return original_hash_band(*args)

    monkeypatch.setattr(kmer_pipeline, "_hash_fasta_band", counted_hash_band)
    first_plan = SequenceBandPlan(len(sequence), 4, 10, (6,))
    list(
        iter_hashed_fasta_bands(
            str(fasta_path),
            "chr1",
            first_plan,
            use_process=False,
            cache_dir=str(cache_dir),
        )
    )
    first_run_calls = calls
    assert first_run_calls == first_plan.band_count
    cached_hashes = load_cached_sequence_hashes(
        str(cache_dir), str(fasta_path), "chr1", len(sequence), 4
    )
    assert isinstance(cached_hashes, np.memmap)
    assert np.array_equal(cached_hashes, expected)

    second_plan = SequenceBandPlan(len(sequence), 4, 7, (5, 9))
    second_bands = list(
        iter_hashed_fasta_bands(
            str(fasta_path),
            "chr1",
            second_plan,
            use_process=False,
            cache_dir=str(cache_dir),
        )
    )

    assert calls == first_run_calls
    assert len(list(cache_dir.glob("*.npy"))) == 1
    for request, band in zip(second_plan.requests(), second_bands):
        assert np.array_equal(
            band.hashes,
            expected[request.start : request.start + request.hash_count],
        )


def test_cache_finalization_failure_does_not_duplicate_bands(tmp_path, monkeypatch):
    fasta_path, sequence = _indexed_fasta(tmp_path)
    plan = SequenceBandPlan(len(sequence), 4, 10, (6,))
    monkeypatch.setattr(
        kmer_pipeline.os,
        "replace",
        lambda *args: (_ for _ in ()).throw(OSError("cache unavailable")),
    )

    bands = list(
        iter_hashed_fasta_bands(
            str(fasta_path),
            "chr1",
            plan,
            use_process=False,
            cache_dir=str(tmp_path / "cache"),
        )
    )

    assert [band.index for band in bands] == list(range(plan.band_count))
