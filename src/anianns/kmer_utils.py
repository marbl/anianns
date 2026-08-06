import mmh3
import pysam
from alive_progress import alive_it
import sys
from typing import Iterable, List, Sequence
import numpy as np
from numba import njit, types
from numba.typed import Dict as NumbaDict

tab_b = bytes.maketrans(b"ACTG", b"TGAC")


@njit(cache=True)
def calculate_hash_distances(hashes):
    """Return distances between consecutive occurrences of each hash."""
    last_indices = NumbaDict.empty(
        key_type=types.int32,
        value_type=types.int64,
    )
    distances = np.empty(len(hashes), dtype=np.int64)
    distance_count = 0
    for index in range(len(hashes)):
        value = hashes[index]
        if value in last_indices:
            distances[distance_count] = index - last_indices[value]
            distance_count += 1
        last_indices[value] = index
    return distances[:distance_count].copy()


def _progress_settings(n: int, k: int):
    total_kmers = n - k + 1
    if total_kmers <= 0:
        return 0, 1
    return total_kmers, max(1, round(n / 77))


def remove_ambiguous_bases(mod_list, k):
    # Ambiguous IUPAC codes
    bases_to_remove = ["R", "Y", "M", "K", "S", "W", "H", "B", "V", "D", "N"]
    kmers_to_remove = set()
    for i in range(len(bases_to_remove)):
        result_string = str(bases_to_remove[i]) * k
        kmers_to_remove.add(mmh3.hash(result_string))
    mod_set = set(mod_list)
    # Remove homopolymers of ambiguous nucleotides
    mod_set.difference_update(kmers_to_remove)
    return mod_set


def convert_set_list_to_sorted_arrays(set_list):
    return [np.array(sorted(s), dtype=np.int32) for s in set_list]


@njit(cache=True)
def _selected_unique(values, selected, start, end):
    """Return sorted unique selected values from one positional slice."""
    start = max(0, start)
    end = min(end, len(values))
    count = 0
    for index in range(start, end):
        if selected[index]:
            count += 1

    result = np.empty(count, dtype=np.int32)
    output_index = 0
    for index in range(start, end):
        if selected[index]:
            result[output_index] = values[index]
            output_index += 1

    if count <= 1:
        return result

    result.sort()
    unique_count = 1
    for index in range(1, count):
        if result[index] != result[unique_count - 1]:
            result[unique_count] = result[index]
            unique_count += 1
    return result[:unique_count].copy()


@njit(cache=True)
def _build_kmer_sets(values, selected, max_len, window, interval):
    non_sets = []
    overlap_sets = []
    for index in range(max_len - 1):
        start = index * window
        end = start + window
        overlap_start = max(0, start - interval)
        overlap_end = end + interval
        non_sets.append(_selected_unique(values, selected, start, end))
        overlap_sets.append(
            _selected_unique(values, selected, overlap_start, overlap_end)
        )
    return overlap_sets, non_sets


def build_kmer_sets(kmer_list, max_len, window, interval, prepend=None, sketch=4):
    """Build sorted window sketches after computing the modulo mask once."""
    if sketch not in (2, 4):
        raise ValueError("sketch must be either 2 or 4")

    values = np.asarray(kmer_list, dtype=np.int32)
    selected = (values != 0) & (values % sketch == 0)
    overlap_sets, non_sets = _build_kmer_sets(
        values, selected, max_len, window, interval
    )

    if prepend is not None and non_sets:
        prepend_values = np.asarray(prepend, dtype=np.int32)
        prepend_values = prepend_values[
            (prepend_values != 0) & (prepend_values % sketch == 0)
        ]
        non_sets[0] = np.unique(np.concatenate((prepend_values, non_sets[0]))).astype(
            np.int32, copy=False
        )

    return overlap_sets, non_sets


def build_kmer_sets_multi(kmer_list, window_configs, sketch=4):
    """Build several window resolutions while calculating the sketch mask once.

    ``window_configs`` maps each window size to
    ``(hash_count, max_len, interval)``. All views must be prefixes of the
    supplied shared band hashes.
    """
    if sketch not in (2, 4):
        raise ValueError("sketch must be either 2 or 4")

    values = np.asarray(kmer_list, dtype=np.int32)
    selected = (values != 0) & (values % sketch == 0)
    results = {}
    for window, (hash_count, max_len, interval) in window_configs.items():
        hash_count = int(hash_count)
        if hash_count < 0 or hash_count > len(values):
            raise ValueError(
                f"window {window} requested {hash_count} hashes from a "
                f"{len(values)}-hash band"
            )
        results[int(window)] = _build_kmer_sets(
            values[:hash_count],
            selected[:hash_count],
            int(max_len),
            int(window),
            int(interval),
        )
    return results


def read_sequence_kmers_from_file(
    filename: str, seqid: str, ksize: int, quiet: bool
) -> List[List[int]]:
    """
    Given a filename and an integer k, returns a list of all k-mers found in the sequences in the file.
    """
    all_kmers = []
    seq = pysam.FastaFile(filename)

    print(f"Retrieving k-mers from {seqid}.... \n")
    kmers_for_seq = []
    for kmer_hash in generate_kmers_from_fasta(seq.fetch(seqid), ksize, quiet):
        kmers_for_seq.append(kmer_hash)
    all_kmers.append(kmers_for_seq)
    print(f"\n{seqid} k-mers retrieved! \n")

    return all_kmers


def generate_kmers_from_fasta(seq: Sequence[str], k: int, quiet: bool) -> Iterable[int]:
    n = len(seq)
    total_kmers, _progress_thresholds = _progress_settings(n, k)
    if total_kmers <= 0:
        return
    indices = range(total_kmers)
    if not quiet:
        indices = alive_it(
            indices,
            title="Hashing k-mers",
            unit=" k-mer",
            enrich_print=False,
            file=sys.stdout,
        )

    bases_to_remove = ["R", "Y", "M", "K", "S", "W", "H", "B", "V", "D", "N"]
    for i in indices:
        # Remove case sensitivity
        kmer = seq[i : i + k].upper()
        # Skip kmer if it contains any ambiguous base
        if any(base in kmer for base in bases_to_remove):
            yield 0

        else:
            fh = mmh3.hash(kmer, seed=42)

            # Calculate reverse complement hash directly without the need for translation
            rc = mmh3.hash(kmer[::-1].translate(tab_b), seed=42)

            yield fh if fh < rc else rc


def generate_kmers_from_fasta_forward_only(
    seq: Sequence[str], k: int, quiet: bool
) -> Iterable[int]:
    n = len(seq)
    total_kmers, _progress_thresholds = _progress_settings(n, k)
    if total_kmers <= 0:
        return
    indices = range(total_kmers)
    if not quiet:
        indices = alive_it(
            indices,
            title="Hashing forward k-mers",
            unit=" k-mer",
            enrich_print=False,
            file=sys.stdout,
        )

    for i in indices:
        # Remove case sensitivity
        kmer = seq[i : i + k].upper()
        fh = mmh3.hash(kmer, seed=42)

        yield fh


def generate_kmers_from_fasta_reverse_only(
    seq: Sequence[str], k: int, quiet: bool
) -> Iterable[int]:
    n = len(seq)
    total_kmers, _progress_thresholds = _progress_settings(n, k)
    if total_kmers <= 0:
        return
    indices = range(total_kmers)
    if not quiet:
        indices = alive_it(
            indices,
            title="Hashing reverse k-mers",
            unit=" k-mer",
            enrich_print=False,
            file=sys.stdout,
        )

    for i in indices:
        # Remove case sensitivity
        kmer = seq[i : i + k].upper()
        rc = mmh3.hash(kmer[::-1].translate(tab_b), seed=42)

        yield rc
