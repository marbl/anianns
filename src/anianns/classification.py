import gzip
import struct
import os
from collections import Counter
from typing import Dict, Tuple, Set, List, Optional


def load_kmer_sets_shared_k(path: str) -> tuple[int, dict[str, set[int]]]:
    with gzip.open(path, "rb") as f:
        k = struct.unpack("<B", f.read(1))[0]
        num_sets = struct.unpack("<H", f.read(2))[0]

        result = {}
        for _ in range(num_sets):
            name_len = struct.unpack("<H", f.read(2))[0]
            name = f.read(name_len).decode("utf-8")
            num_kmers = struct.unpack("<I", f.read(4))[0]
            kmers = {struct.unpack("<i", f.read(4))[0] for _ in range(num_kmers)}
            result[name] = kmers

    return k, result


def load_all_kmer_dbs(directory: str) -> Dict[str, Tuple[int, Dict[str, Set[int]]]]:
    """
    Load all .db files from a directory using load_kmer_sets_shared_k.

    Args:
        directory: Path to directory containing .db files

    Returns:
        Dictionary mapping filename (without .db extension) to (k, kmer_sets) tuples

    Raises:
        FileNotFoundError: If directory doesn't exist
        ValueError: If no .db files found in directory
    """
    if not os.path.exists(directory):
        raise FileNotFoundError(f"Directory not found: {directory}")

    if not os.path.isdir(directory):
        raise ValueError(f"Path is not a directory: {directory}")

    db_files = [f for f in os.listdir(directory) if f.endswith(".db")]

    if not db_files:
        raise ValueError(f"No .db files found in directory: {directory}")

    result = {}
    for db_file in db_files:
        file_path = os.path.join(directory, db_file)
        # Use filename without .db extension as key
        key = os.path.splitext(db_file)[0]
        try:
            k, kmer_sets = load_kmer_sets_shared_k(file_path)
            result[key] = (k, kmer_sets)
        except Exception as e:
            print(f"Warning: Failed to load {db_file}: {e}")
            continue

    return result


def classify_kmers(
    query_kmers: Set[int],
    kmer_supersets: Dict[str, Set[int]],
    min_overlap: int = 1,
    verbose: bool = True,
) -> Tuple[Optional[str], List[Tuple[str, int, float]]]:
    """
    Classify a set of k-mers against known k-mer databases.

    Args:
        query_kmers: Set of k-mers to classify
        kmer_supersets: Dictionary of {database_name: combined_kmer_set}
        min_overlap: Minimum overlap required for classification (default: 1)
        verbose: Print detailed results (default: False)

    Returns:
        Tuple of:
        - Best match database name (None if no match meets min_overlap)
        - List of (database_name, overlap_count, overlap_percentage) sorted by overlap_count descending
    """
    if not query_kmers:
        return "Unknown", []

    results = []
    query_size = len(query_kmers)

    for db_name, db_kmers in kmer_supersets.items():
        overlap = query_kmers & db_kmers
        overlap_count = len(overlap)
        overlap_percentage = (
            (overlap_count / query_size) * 100 if query_size > 0 else 0.0
        )

        results.append((db_name, overlap_count, overlap_percentage))

    # Sort by overlap count (descending)
    results.sort(key=lambda x: x[1], reverse=True)

    # Show top 3 results if verbose
    verbose = False
    if verbose:
        print("Top 3 classification matches:")
        for i, (db_name, overlap_count, overlap_percentage) in enumerate(results[:3]):
            print(
                f"  {i+1}. {db_name}: {overlap_count} overlap ({overlap_percentage:.2f}%)"
            )

    best_match = "Unknown"
    if results:
        top_hit = results[0]
        if top_hit[1] >= min_overlap and top_hit[2] > 25.0:
            best_match = top_hit[0]

    return best_match, results


def calculate_distances(numbers):
    indices_map = {}
    distances = []

    # Iterate through the list and store indices of each number
    for idx, num in enumerate(numbers):
        if num in indices_map:
            # Calculate the distance between current and previous occurrence
            prev_idx = indices_map[num][-1]
            distance = idx - prev_idx
            distances.append(distance)
            # Update the list of indices for this number
            indices_map[num].append(idx)
        else:
            # If the number is seen for the first time, just store its index
            indices_map[num] = [idx]

    return distances


def top_n_frequent_distances(distances, n=5):
    # Count the occurrences of each distance
    distance_counts = Counter(distances)
    # Get the top n most common distances
    top_n = distance_counts.most_common(n)

    # Print the results
    print(f"Top {n} distances and their counts:")
    for distance, count in top_n:
        print(f"Distance: {distance}, Count: {count}")
