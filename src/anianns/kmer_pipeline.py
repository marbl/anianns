"""Streaming FASTA band planning and canonical k-mer production."""

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
import hashlib
import math
import multiprocessing
import os
from pathlib import Path
import tempfile
from typing import Dict, Iterable, Iterator, Sequence

import numpy as np
import pysam
from numba import njit, prange, set_num_threads

from anianns.kmer_utils import tab_b


HASH_CACHE_VERSION = 1


@dataclass(frozen=True)
class WindowPlan:
    """Derived settings for one matrix window size."""

    window: int
    interval: int
    max_len: int


@dataclass(frozen=True)
class BandRequest:
    """A zero-based k-mer interval to fetch and hash."""

    index: int
    start: int
    hash_count: int


@dataclass(frozen=True)
class HashedBand:
    """Canonical hashes returned by a band producer."""

    index: int
    start: int
    hashes: Sequence[int]


class SequenceBandPlan:
    """
    Plan shared k-mer batches for one or more prospective window sizes.

    Each FASTA band is fetched and hashed once with enough right-hand overlap
    for the largest requested window, and each resolution consumes only its
    own view of that shared hash stream.
    """

    def __init__(
        self,
        sequence_length: int,
        kmer: int,
        band_height: int,
        windows: Sequence[int],
    ):
        if sequence_length < 0:
            raise ValueError("sequence_length must be non-negative")
        if kmer <= 0:
            raise ValueError("kmer must be positive")
        if band_height <= 0:
            raise ValueError("band_height must be positive")
        if not windows or any(window <= 0 for window in windows):
            raise ValueError("at least one positive window size is required")

        self.sequence_length = sequence_length
        self.kmer = kmer
        self.band_height = band_height
        self.windows = tuple(dict.fromkeys(windows))
        self.window_plans: Dict[int, WindowPlan] = {
            window: WindowPlan(
                window=window,
                interval=(window + 1) // 2,
                max_len=(band_height + window) // window,
            )
            for window in self.windows
        }
        self.total_kmers = max(0, sequence_length - kmer + 1)
        self.max_interval = max(plan.interval for plan in self.window_plans.values())
        self.band_count = (
            math.ceil(self.total_kmers / band_height) if self.total_kmers else 0
        )

    def requests(self) -> Iterator[BandRequest]:
        """Yield complete, non-growing band requests through sequence end."""
        for index in range(self.band_count):
            start = index * self.band_height
            available = self.total_kmers - start
            hash_count = min(self.band_height + self.max_interval, available)
            yield BandRequest(index=index, start=start, hash_count=hash_count)

    def hashes_for_window(self, band: HashedBand, window: int) -> Sequence[int]:
        """Return the shared batch prefix required by one window size."""
        try:
            window_plan = self.window_plans[window]
        except KeyError as error:
            raise ValueError(
                f"window size {window} is not part of this plan"
            ) from error

        available = self.total_kmers - band.start
        required = min(self.band_height + window_plan.interval, available)
        if len(band.hashes) < required:
            raise ValueError(
                f"band {band.index} contains {len(band.hashes)} hashes; "
                f"window {window} requires {required}"
            )
        return band.hashes[:required]


@njit(cache=True, inline="always")
def _rotate_left_32(value, shift):
    return np.uint32((value << shift) | (value >> (32 - shift)))


@njit(cache=True, inline="always")
def _finalize_murmurhash3_32(value):
    value ^= value >> 16
    value = np.uint32(value * np.uint32(0x85EBCA6B))
    value ^= value >> 13
    value = np.uint32(value * np.uint32(0xC2B2AE35))
    value ^= value >> 16
    return value


@njit(cache=True, parallel=True)
def _forward_murmurhash3_windows(sequence_bytes, kmer, seed):
    total_kmers = len(sequence_bytes) - kmer + 1
    hashes = np.empty(max(0, total_kmers), dtype=np.int32)
    c1 = np.uint32(0xCC9E2D51)
    c2 = np.uint32(0x1B873593)

    for start in prange(total_kmers):
        value = np.uint32(seed)
        block_count = kmer // 4
        for block in range(block_count):
            offset = start + (block * 4)
            block_value = np.uint32(
                np.uint32(sequence_bytes[offset])
                | (np.uint32(sequence_bytes[offset + 1]) << 8)
                | (np.uint32(sequence_bytes[offset + 2]) << 16)
                | (np.uint32(sequence_bytes[offset + 3]) << 24)
            )
            block_value = np.uint32(block_value * c1)
            block_value = _rotate_left_32(block_value, 15)
            block_value = np.uint32(block_value * c2)
            value ^= block_value
            value = _rotate_left_32(value, 13)
            value = np.uint32(value * np.uint32(5) + np.uint32(0xE6546B64))

        tail_offset = start + (block_count * 4)
        tail_size = kmer & 3
        tail = np.uint32(0)
        if tail_size == 3:
            tail ^= np.uint32(sequence_bytes[tail_offset + 2]) << 16
        if tail_size >= 2:
            tail ^= np.uint32(sequence_bytes[tail_offset + 1]) << 8
        if tail_size >= 1:
            tail ^= np.uint32(sequence_bytes[tail_offset])
            tail = np.uint32(tail * c1)
            tail = _rotate_left_32(tail, 15)
            tail = np.uint32(tail * c2)
            value ^= tail

        value ^= np.uint32(kmer)
        value = _finalize_murmurhash3_32(value)
        if value >= np.uint32(0x80000000):
            hashes[start] = np.int32(np.int64(value) - 0x100000000)
        else:
            hashes[start] = np.int32(value)

    return hashes


@njit(cache=True, inline="always")
def _is_ambiguous_base(value):
    """Match the historical IUPAC ambiguity filter exactly."""
    return (
        value == ord("R")
        or value == ord("Y")
        or value == ord("M")
        or value == ord("K")
        or value == ord("S")
        or value == ord("W")
        or value == ord("H")
        or value == ord("B")
        or value == ord("V")
        or value == ord("D")
        or value == ord("N")
    )


@njit(cache=True)
def _canonicalize_forward_hashes(
    sequence_bytes, forward_hashes, reverse_complement_hashes, kmer
):
    """Canonicalize two forward-hash vectors and mask ambiguous windows."""
    total_kmers = len(forward_hashes)
    ambiguous_count = 0
    for index in range(kmer):
        if _is_ambiguous_base(sequence_bytes[index]):
            ambiguous_count += 1

    for index in range(total_kmers):
        if ambiguous_count:
            forward_hashes[index] = 0
        else:
            reverse_hash = reverse_complement_hashes[total_kmers - index - 1]
            if reverse_hash < forward_hashes[index]:
                forward_hashes[index] = reverse_hash

        if index + 1 < total_kmers:
            if _is_ambiguous_base(sequence_bytes[index]):
                ambiguous_count -= 1
            if _is_ambiguous_base(sequence_bytes[index + kmer]):
                ambiguous_count += 1
    return forward_hashes


def canonical_kmer_hashes(sequence: str, kmer: int) -> np.ndarray:
    """
    Batch canonical mmh3 hashing while preserving existing hash semantics.

    Every forward k-mer and reverse-complement k-mer is hashed in two parallel
    Numba passes. Reversing the reverse-complement hash vector aligns it with
    the original sequence, after which a linear pass selects the canonical
    minimum and applies the historical ambiguous-base mask.
    """
    if kmer <= 0:
        raise ValueError("kmer must be positive")

    sequence = sequence.upper()
    total_kmers = len(sequence) - kmer + 1
    if total_kmers <= 0:
        return np.empty(0, dtype=np.int32)

    sequence_bytes = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
    reverse_complement = sequence[::-1].translate(tab_b)
    reverse_complement_bytes = np.frombuffer(
        reverse_complement.encode("ascii"), dtype=np.uint8
    )
    forward_hashes = _forward_murmurhash3_windows(sequence_bytes, kmer, 42)
    reverse_complement_hashes = _forward_murmurhash3_windows(
        reverse_complement_bytes, kmer, 42
    )
    return _canonicalize_forward_hashes(
        sequence_bytes,
        forward_hashes,
        reverse_complement_hashes,
        kmer,
    )


def forward_kmer_hashes(sequence: str, kmer: int) -> np.ndarray:
    """Batch forward-only mmh3 hashes without allocating k-mer substrings."""
    if kmer <= 0:
        raise ValueError("kmer must be positive")
    encoded = np.frombuffer(sequence.upper().encode("ascii"), dtype=np.uint8)
    if len(encoded) < kmer:
        return np.empty(0, dtype=np.int32)
    return _forward_murmurhash3_windows(encoded, kmer, 42)


_FASTA_HANDLE_CACHE = {}


def _hash_fasta_band(
    fasta_path: str,
    seq_id: str,
    sequence_length: int,
    kmer: int,
    request: BandRequest,
    worker_threads: int = None,
) -> HashedBand:
    """Fetch and hash one band. This top-level function is process-picklable."""
    if worker_threads is not None:
        set_num_threads(worker_threads)
    fasta = _FASTA_HANDLE_CACHE.get(fasta_path)
    if fasta is None:
        fasta = pysam.FastaFile(fasta_path)
        _FASTA_HANDLE_CACHE[fasta_path] = fasta

    base_end = min(sequence_length, request.start + request.hash_count + kmer - 1)
    sequence = fasta.fetch(seq_id, request.start, base_end)
    hashes = canonical_kmer_hashes(sequence, kmer)
    if len(hashes) != request.hash_count:
        raise ValueError(
            f"Expected {request.hash_count} hashes for band {request.index}, "
            f"received {len(hashes)}"
        )
    return HashedBand(index=request.index, start=request.start, hashes=hashes)


def hash_cache_path(
    cache_dir: str, fasta_path: str, seq_id: str, sequence_length: int, kmer: int
) -> Path:
    """Return the content-versioned cache path for one sequence and k-mer size."""
    stat = os.stat(fasta_path)
    identity = "\0".join(
        (
            str(HASH_CACHE_VERSION),
            str(Path(fasta_path).resolve()),
            str(stat.st_size),
            str(stat.st_mtime_ns),
            seq_id,
            str(sequence_length),
            str(kmer),
        )
    )
    digest = hashlib.sha256(identity.encode("utf-8")).hexdigest()
    return Path(cache_dir) / f"{digest}.npy"


def _load_cached_hashes(cache_path: Path, total_kmers: int):
    try:
        hashes = np.load(cache_path, mmap_mode="r", allow_pickle=False)
    except (OSError, ValueError):
        return None
    if hashes.dtype != np.int32 or hashes.shape != (total_kmers,):
        return None
    return hashes


def load_cached_sequence_hashes(
    cache_dir: str,
    fasta_path: str,
    seq_id: str,
    sequence_length: int,
    kmer: int,
):
    """Load the complete canonical hash vector when a valid cache is present."""
    if cache_dir is None:
        return None
    try:
        cache_path = hash_cache_path(
            cache_dir, fasta_path, seq_id, sequence_length, kmer
        )
    except OSError:
        return None
    return _load_cached_hashes(cache_path, max(0, sequence_length - kmer + 1))


def _iter_uncached_hashed_fasta_bands(
    fasta_path: str,
    seq_id: str,
    plan: SequenceBandPlan,
    *,
    use_process: bool,
    worker_threads: int,
) -> Iterable[HashedBand]:
    """Produce hashed bands without consulting or writing the disk cache."""
    requests = iter(plan.requests())
    try:
        first_request = next(requests)
    except StopIteration:
        return

    worker_args = (fasta_path, seq_id, plan.sequence_length, plan.kmer)
    if not use_process or plan.band_count == 1:
        yield _hash_fasta_band(*worker_args, first_request)
        for request in requests:
            yield _hash_fasta_band(*worker_args, request)
        return

    context = multiprocessing.get_context("spawn")
    try:
        executor = ProcessPoolExecutor(max_workers=1, mp_context=context)
    except (OSError, NotImplementedError):
        # Restricted runtimes may not expose the semaphore/sysconf facilities
        # required by ProcessPoolExecutor. Preserve correctness by streaming
        # synchronously rather than failing annotation startup.
        yield _hash_fasta_band(*worker_args, first_request)
        for request in requests:
            yield _hash_fasta_band(*worker_args, request)
        return
    try:
        current_future = executor.submit(
            _hash_fasta_band, *worker_args, first_request, worker_threads
        )
        while current_future is not None:
            band = current_future.result()
            try:
                next_request = next(requests)
            except StopIteration:
                next_future = None
            else:
                # Submit before yielding so hashing overlaps parent matrix work.
                next_future = executor.submit(
                    _hash_fasta_band, *worker_args, next_request, worker_threads
                )
            yield band
            current_future = next_future
    finally:
        executor.shutdown(wait=True, cancel_futures=True)


def iter_hashed_fasta_bands(
    fasta_path: str,
    seq_id: str,
    plan: SequenceBandPlan,
    *,
    use_process: bool = True,
    cache_dir: str = None,
    worker_threads: int = 1,
) -> Iterable[HashedBand]:
    """
    Yield bands in order with at most one future band queued ahead.

    A single producer process avoids the GIL during canonical hashing while the
    parent performs Numba matrix work. The one-band bound caps additional hash
    memory and avoids oversubscribing Numba's parallel matrix workers.
    """
    if plan.total_kmers == 0:
        return

    if cache_dir is None:
        yield from _iter_uncached_hashed_fasta_bands(
            fasta_path,
            seq_id,
            plan,
            use_process=use_process,
            worker_threads=worker_threads,
        )
        return

    try:
        cache_path = hash_cache_path(
            cache_dir, fasta_path, seq_id, plan.sequence_length, plan.kmer
        )
        cache_path.parent.mkdir(parents=True, exist_ok=True)
    except OSError:
        yield from _iter_uncached_hashed_fasta_bands(
            fasta_path,
            seq_id,
            plan,
            use_process=use_process,
            worker_threads=worker_threads,
        )
        return

    cached_hashes = _load_cached_hashes(cache_path, plan.total_kmers)
    if cached_hashes is not None:
        for request in plan.requests():
            yield HashedBand(
                index=request.index,
                start=request.start,
                hashes=cached_hashes[
                    request.start : request.start + request.hash_count
                ],
            )
        return

    try:
        temporary = tempfile.NamedTemporaryFile(
            dir=cache_path.parent,
            prefix=f".{cache_path.stem}.",
            suffix=".npy",
            delete=False,
        )
    except OSError:
        yield from _iter_uncached_hashed_fasta_bands(
            fasta_path,
            seq_id,
            plan,
            use_process=use_process,
            worker_threads=worker_threads,
        )
        return
    temporary_path = Path(temporary.name)
    temporary.close()
    try:
        cache_hashes = np.lib.format.open_memmap(
            temporary_path,
            mode="w+",
            dtype=np.int32,
            shape=(plan.total_kmers,),
        )
    except (OSError, ValueError):
        try:
            temporary_path.unlink()
        except OSError:
            pass
        yield from _iter_uncached_hashed_fasta_bands(
            fasta_path,
            seq_id,
            plan,
            use_process=use_process,
            worker_threads=worker_threads,
        )
        return

    cache_writable = True
    try:
        for band in _iter_uncached_hashed_fasta_bands(
            fasta_path,
            seq_id,
            plan,
            use_process=use_process,
            worker_threads=worker_threads,
        ):
            core_count = min(plan.band_height, plan.total_kmers - band.start)
            try:
                if cache_writable:
                    cache_hashes[band.start : band.start + core_count] = band.hashes[
                        :core_count
                    ]
            except (OSError, ValueError):
                cache_writable = False
            yield band

        if cache_writable:
            try:
                cache_hashes.flush()
                del cache_hashes
                cache_hashes = None
                os.replace(temporary_path, cache_path)
            except (OSError, ValueError):
                # Caching is an optimization and must never fail annotation.
                pass
    finally:
        if cache_hashes is not None:
            try:
                cache_hashes.flush()
            except (OSError, ValueError):
                pass
            del cache_hashes
        if temporary_path.exists():
            try:
                temporary_path.unlink()
            except OSError:
                pass
