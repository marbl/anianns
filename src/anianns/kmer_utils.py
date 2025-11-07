import mmh3
import pysam
from typing import Iterable, List, Sequence
import numpy as np

tab_b = bytes.maketrans(b"ACTG", b"TGAC")

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

def build_kmer_sets(kmer_list, max_len, window, interval, prepend=None):
    non_sets = []
    overlap_sets = []

    for i in range(max_len - 1):
        start = i * window
        end = start + window
        ostart = max(0, start - interval)
        oend = end + interval

        # non-overlap slice (prepend only on the very first window if requested)
        if prepend is not None and ostart == 0:
            seq_non = prepend + kmer_list[start:end]
        else:
            seq_non = kmer_list[start:end]

        # build your sets in one pass each, remove 0 k-mers
        non_sets.append({x for x in seq_non if x % 4 == 0 and x != 0})
        overlap_sets.append({x for x in kmer_list[ostart:oend] if x % 4 == 0 and x != 0})

    # return exactly as before (overlap first, then non-overlap)
    return (
        convert_set_list_to_sorted_arrays(overlap_sets),
        convert_set_list_to_sorted_arrays(non_sets),
    )

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
    if not quiet:
        progress_thresholds = round(n / 77)
        print_progress_bar(
            0, n - k + 1, prefix="Progress:", suffix="Complete", length=40
        )

    bases_to_remove = ["R", "Y", "M", "K", "S", "W", "H", "B", "V", "D", "N"]
    for i in range(n - k + 1):
        if not quiet:
            if i % progress_thresholds == 0:
                print_progress_bar(
                    i, n - k + 1, prefix="Progress:", suffix="Complete", length=40
                )
            if i == n - k:
                print_progress_bar(
                    n - k + 1,
                    n - k + 1,
                    prefix="Progress:",
                    suffix="Completed",
                    length=40,
                )
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
    if not quiet:
        progress_thresholds = round(n / 77)
        print_progress_bar(
            0, n - k + 1, prefix="Progress:", suffix="Complete", length=40
        )

    for i in range(n - k + 1):
        if not quiet:
            if i % progress_thresholds == 0:
                print_progress_bar(
                    i, n - k + 1, prefix="Progress:", suffix="Complete", length=40
                )
            if i == n - k:
                print_progress_bar(
                    n - k + 1,
                    n - k + 1,
                    prefix="Progress:",
                    suffix="Completed",
                    length=40,
                )
        # Remove case sensitivity
        kmer = seq[i : i + k].upper()
        fh = kmer

        yield fh


def generate_kmers_from_fasta_reverse_only(
    seq: Sequence[str], k: int, quiet: bool
) -> Iterable[int]:
    n = len(seq)
    if not quiet:
        progress_thresholds = round(n / 77)
        print_progress_bar(
            0, n - k + 1, prefix="Progress:", suffix="Complete", length=40
        )

    for i in range(n - k + 1):
        if not quiet:
            if i % progress_thresholds == 0:
                print_progress_bar(
                    i, n - k + 1, prefix="Progress:", suffix="Complete", length=40
                )
            if i == n - k:
                print_progress_bar(
                    n - k + 1,
                    n - k + 1,
                    prefix="Progress:",
                    suffix="Completed",
                    length=40,
                )
        # Remove case sensitivity
        kmer = seq[i : i + k].upper()
        rc = mmh3.hash(kmer[::-1].translate(tab_b), seed=42)

        yield rc


def print_progress_bar(
    iteration,
    total,
    prefix="",
    suffix="",
    decimals=1,
    length=100,
    fill="█",
    printEnd="\r",
):
    percent = f"{100 * (iteration / total):.{decimals}f}"
    filledLength = int(length * iteration // total)
    bar = [fill] * filledLength + ["-"] * (length - filledLength)
    bar_str = "".join(bar)
    print(f"\r{prefix} |{bar_str}| {percent}% {suffix}", end=printEnd)
    if iteration == total:
        print()
