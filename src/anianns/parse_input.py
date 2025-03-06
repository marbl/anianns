from typing import Iterable, List, Sequence
from moddotplot.parse_fasta import generateKmersFromFasta
import pysam
import sys
import mmh3

tab_b = bytes.maketrans(b"ACTG", b"TGAC")


def printProgressBar(
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


def readKmersFromFileNonHashed(
    filename: str, ksize: int, quiet: bool
) -> List[List[int]]:
    """
    Given a filename and an integer k, returns a list of all k-mers found in the sequences in the file.
    """
    all_kmers = []
    seq = pysam.FastaFile(filename)

    for seq_id in seq.references:
        print(f"Retrieving k-mers from {seq_id}.... \n")
        kmers_for_seq = []
        for kmer_hash in generateKmersFromFastaNonHashed(
            seq.fetch(seq_id), ksize, quiet
        ):
            kmers_for_seq.append(kmer_hash)
        all_kmers.append(kmers_for_seq)
        print(f"\n{seq_id} k-mers retrieved! \n")

    return all_kmers


def generateKmersFromFastaNonHashed(
    seq: Sequence[str], k: int, quiet: bool
) -> Iterable[int]:
    n = len(seq)
    if not quiet:
        progress_thresholds = round(n / 77)
        printProgressBar(0, n - k + 1, prefix="Progress:", suffix="Complete", length=40)

    for i in range(n - k + 1):
        if not quiet:
            if i % progress_thresholds == 0:
                printProgressBar(
                    i, n - k + 1, prefix="Progress:", suffix="Complete", length=40
                )
            if i == n - k:
                printProgressBar(
                    n - k + 1,
                    n - k + 1,
                    prefix="Progress:",
                    suffix="Completed",
                    length=40,
                )
        # Remove case sensitivity
        kmer = seq[i : i + k].upper()

        yield kmer


def readSequenceKmersFromFile(
    filename: str, seqid: str, ksize: int, quiet: bool
) -> List[List[int]]:
    """
    Given a filename and an integer k, returns a list of all k-mers found in the sequences in the file.
    """
    all_kmers = []
    seq = pysam.FastaFile(filename)

    print(f"Retrieving k-mers from {seqid}.... \n")
    kmers_for_seq = []
    for kmer_hash in generateKmersFromFasta(seq.fetch(seqid), ksize, quiet):
        kmers_for_seq.append(kmer_hash)
    all_kmers.append(kmers_for_seq)
    print(f"\n{seqid} k-mers retrieved! \n")

    return all_kmers
