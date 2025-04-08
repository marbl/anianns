from typing import Iterable, List, Sequence, Callable
from moddotplot.parse_fasta import generateKmersFromFasta
import pysam
import sys
import mmh3
import struct
import zlib
from anianns.const import HEADER_FORMAT, HEADER_SIZE 
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

def readKmersFromFileCustom(
    filename: str, ksize: int, quiet: bool, kmer_fn: Callable[[str, int, bool], List[int]]
) -> List[List[int]]:
    """
    Given a filename and an integer k, returns a list of all k-mers found in the sequences in the file.
    """
    all_kmers = []
    seq = pysam.FastaFile(filename)

    for seq_id in seq.references:
        print(f"Retrieving k-mers from {seq_id}.... \n")
        kmers_for_seq = []
        for kmer_hash in kmer_fn(
            seq.fetch(seq_id), ksize, quiet
        ):
            kmers_for_seq.append(kmer_hash)
        all_kmers.append(kmers_for_seq)
        print(f"\n{seq_id} k-mers retrieved! \n")

    return all_kmers

def readKmersFromSeqCustom(
    filename, ksize: int, quiet: bool, kmer_fn: Callable[[str, int, bool], List[int]]
) -> List[List[int]]:
    """
    Given a filename and an integer k, returns a list of all k-mers found in the sequences in the file.
    """
    kmers_for_seq = []
    for kmer_hash in kmer_fn(
        filename, ksize, quiet
    ):
        kmers_for_seq.append(kmer_hash)

    return kmers_for_seq


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

def generateKmersFromFastaForward(seq: Sequence[str], k: int, quiet: bool) -> Iterable[int]:
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
        fh = mmh3.hash(kmer, seed=42)

        yield fh

def generateKmersFromFastaReverse(seq: Sequence[str], k: int, quiet: bool) -> Iterable[int]:
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

        # Calculate reverse complement hash directly without the need for translation
        rc = mmh3.hash(kmer[::-1].translate(tab_b), seed=42)

        yield rc


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

def readSequenceKmersFromFileForward(
    filename: str, seqid: str, ksize: int, quiet: bool
) -> List[List[int]]:
    """
    Given a filename and an integer k, returns a list of all k-mers found in the sequences in the file.
    """
    all_kmers = []
    seq = pysam.FastaFile(filename)

    print(f"Retrieving k-mers from {seqid}.... \n")
    kmers_for_seq = []
    for kmer_hash in generateKmersFromFastaForward(seq.fetch(seqid), ksize, quiet):
        kmers_for_seq.append(kmer_hash)
    all_kmers.append(kmers_for_seq)
    print(f"\n{seqid} k-mers retrieved! \n")

    return all_kmers

def readSequenceKmersFromFileReverse(
    filename: str, seqid: str, ksize: int, quiet: bool
) -> List[List[int]]:
    """
    Given a filename and an integer k, returns a list of all k-mers found in the sequences in the file.
    """
    all_kmers = []
    seq = pysam.FastaFile(filename)

    print(f"Retrieving k-mers from {seqid}.... \n")
    kmers_for_seq = []
    for kmer_hash in generateKmersFromFastaReverse(seq.fetch(seqid), ksize, quiet):
        kmers_for_seq.append(kmer_hash)
    all_kmers.append(kmers_for_seq)
    print(f"\n{seqid} k-mers retrieved! \n")

    return all_kmers

# Writing a function that parses a bed file, then extracts the sequence within
def extract_region(fasta_file, chr, region_start, region_end):
    # Open the FASTA file
    fasta = pysam.FastaFile(fasta_file)

    try:
    # Extract sequence from 1-based start to end (pysam is 0-based, so adjust)
        sequence = fasta.fetch(chr, region_start - 1, region_end)
    except:
        print(f"Unable to find {chr}. Make sure name is fasta file = name in index or bed file.\n")
        sequence = None

    # Close FASTA file
    fasta.close()

    return sequence

def read_bedfile(filepath):
    """Efficiently reads a BED-like file and returns parsed data as a list of tuples."""
    bed_data = []
    with open(filepath, "r") as f:
        next(f)
        append = bed_data.append  # Local variable reference for faster access
        for line in f:
            parts = line.split("\t")  # Avoids unnecessary .strip() as split() handles it
            #chrom, coords = parts[0].split(":")
            #start, end = map(int, coords.split("-"))
            append((parts[0], int(parts[1]), int(parts[2]), parts[3], int(parts[4]), parts[5], int(parts[6]), int(parts[7]), parts[8].split("\n")[0]))  # Direct conversion
    return bed_data

def minhash_sketch(hash_list, s=10_000):
    """
    Generate a MinHash sketch of size s from a set of unique hash values.

    Parameters:
    - hash_set (set[int]): A set of unique hash values.
    - s (int): The size of the MinHash sketch (default: 10,000).

    Returns:
    - list[int]: The MinHash sketch containing the s smallest hashes.
    """
    hash_set = set(hash_list)
    sorted_list = sorted(hash_set)[:s]
    return set(sorted_list)

def load_kmers_binary(filename):
    """Load k-mer hashes stored as signed 64-bit integers from a compressed binary file."""
    with open(filename, "rb") as f:
        compressed_data = f.read()

    binary_data = zlib.decompress(compressed_data)
    kmers = [struct.unpack("q", binary_data[i:i+8])[0] for i in range(0, len(binary_data), 8)]

    print(f"Loaded {len(kmers)} k-mers from {filename}")
    return kmers

def save_kmers_binary(kmer_sets, filename, append=False):
    """
    Save multiple k-mer sets into a compressed binary file with headers.

    Parameters:
    - kmer_sets (dict[str, set[int]]): Dictionary mapping set names to k-mer sets.
    - filename (str): File to save the data.
    - append (bool): If True, append to the existing file.
    """
    mode = "ab" if append else "wb"

    with open(filename, mode) as f:
        for set_name, kmers in kmer_sets.items():
            # Ensure name is exactly 32 bytes (padded or truncated)
            name_bytes = set_name.encode("utf-8")[:32].ljust(32, b"\0")

            # Pack the header (name + number of kmers)
            header = struct.pack(HEADER_FORMAT, name_bytes, len(kmers))

            # Pack k-mer data
            kmer_data = b"".join(struct.pack("q", kmer) for kmer in kmers)
            
            # Compress the header + data
            compressed_data = zlib.compress(header + kmer_data)

            # Write to file
            f.write(compressed_data)

            print(f"Saved set '{set_name}' ({len(kmers)} k-mers) to {filename} (compressed size: {len(compressed_data)} bytes, append={append})")

def load_kmers_binary(filename):
    """
    Load multiple named k-mer sets from a compressed binary file.

    Parameters:
    - filename (str): File to load the data from.

    Returns:
    - dict[str, set[int]]: Dictionary mapping set names to their corresponding k-mer sets.
    """
    kmer_sets = {}

    with open(filename, "rb") as f:
        while True:
            # Try reading a compressed block
            try:
                compressed_data = f.read()  # Read all remaining bytes
                if not compressed_data:
                    break

                # Decompress directly (assuming one large compressed file)
                binary_data = zlib.decompress(compressed_data)

                # Ensure header exists
                if len(binary_data) < HEADER_SIZE:
                    print("Error: Incomplete header.")
                    break

                # Extract header
                name_bytes, num_kmers = struct.unpack(HEADER_FORMAT, binary_data[:HEADER_SIZE])
                set_name = name_bytes.rstrip(b"\0").decode("utf-8")  # Strip null padding

                # Read k-mers
                kmers = set()
                for i in range(HEADER_SIZE, len(binary_data), 8):
                    kmers.add(struct.unpack("q", binary_data[i:i+8])[0])

                # Store in dictionary
                kmer_sets[set_name] = kmers
                print(f"Loaded set '{set_name}' ({num_kmers} k-mers)")

                break  # Stop after reading one full compressed block

            except zlib.error as e:
                print(f"Error: Failed to decompress data ({e}).")
                break

    return kmer_sets