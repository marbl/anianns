#from fcntl import F_SEAL_SEAL
from anianns.parse_input import generateKmersFromFastaNonHashed, read_bedfile, extract_region, generateKmersFromFastaForward, generateKmersFromFastaReverse
from itertools import islice
import math
import random
from anianns.parse_input import load_kmers_binary, minhash_sketch
from pathlib import Path
import time
from collections import Counter
from anianns.const import COLOR_MAPPING

def random_multiple(start, end, step):
    first = ((start // step) + 1) * step  # First valid multiple > start
    last = (end // step) * step  # Last valid multiple ≤ end
    
    if first > end:
        raise ValueError("No valid multiples in range.")
    
    return random.choice(range(first, last + 1, step))

# Ideally this should be based on runs, but this works for now
def detect_transitions(data):
    """Identifies points where runs of zeros transition to values > 50 and vice versa."""
    transitions = []
    prev_state = 0 if data[0] == 0 else 1  # 0 for zero-run, 1 for >50-run
    change_potential = False
    change_value = None

    for entry in data[1:]:
        value = entry[1]
        breakpoint = entry[0]
        current_state = 0 if value == 0 else 1  # 0 for zero-run, 1 for >50-run

        if current_state != prev_state and entry != data[1:][0]:
            change_potential = True
            change_value = breakpoint

        elif change_potential:
            transitions.append(change_value)  # Store the range where transition occurs
            change_potential = False
            change_value = None

        prev_state = current_state  # Update previous state

    return transitions

def lists_within_step(list1, list2, step_size):
    if len(list1) != len(list2):
        return F_SEAL_SEAL
    
    pairs = [(min(a, b), max(a, b)) for a, b in zip(list1, list2)]
    
    if all(abs(a - b) <= step_size for a, b in zip(list1, list2)):
        return sorted(pairs)
    return None

def determine_strand(bed_file, step_size, fasta_file, k):
    annotation_entries = read_bedfile(bed_file)

    for entry in annotation_entries:
        chrom, start, end = entry[:3]
        size = end - start

        seq = extract_region(fasta_file, chrom, start, end)
        forward_kmers = generateKmersFromFastaForward(seq, k, True)
        reverse_kmers = generateKmersFromFastaReverse(seq, k, True)

        rv_list = list(islice(reverse_kmers, 1, size))
        query_1 = set(islice(forward_kmers, 1, 2000))
        query_2 = set(islice(forward_kmers, 76001, 78000))

        strand_orientation_1 = []
        strand_orientation_2 = []

        for j in range(math.ceil(size / step_size)):
            window_start = j * step_size + 1
            window_end = min(window_start + step_size, size)

            set_size_q1 = len(query_1.intersection(rv_list[window_start:window_end]))
            set_size_q2 = len(query_2.intersection(rv_list[window_start:window_end]))

            strand_orientation_1.append((window_start, set_size_q1))
            strand_orientation_2.append((window_start, set_size_q2))

        return lists_within_step(
            detect_transitions(strand_orientation_1), 
            detect_transitions(strand_orientation_2), 
            step_size
        )

def find_files_with_suffix(directory: str, suffix: str):
    return list(Path(directory).glob(f"*{suffix}"))

def classify_seq(fasta_file, bed_file, k, sat_directory):
    annotation_entries = read_bedfile(bed_file)

    # Load k-mer database
    sat_db = {}
    print(f"Loading satellite db\n")
    for satellite in sat_directory:
        sat = load_kmers_binary(satellite)
        #sat_db.append(sat)
        for key in sat:
            sat_db[key] = minhash_sketch(sat[key], s=100000)
            #sat_db[key] = sat[key]
    print(f"Satellites loaded! Annotating satellite regions...\n")
    new_annotations = []
    for entry in annotation_entries:
        chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb = entry[:9]
        result = "other"
        #chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb = interval
        size = chromEnd - chromStart
        seq = extract_region(fasta_file, chrom, chromStart, chromEnd)
        forward_kmers = generateKmersFromFastaForward(seq, k, True)
        reverse_kmers = generateKmersFromFastaReverse(seq, k, True)
        query_set = set(list(islice(reverse_kmers, 1, size))) | set(list(islice(forward_kmers, 1, size)))
        lengths = [(key, len(query_set & value)) for key, value in sat_db.items()]
        # Sort by length in descending order
        lengths.sort(key=lambda x: x[1], reverse=True)

        # Check if the largest value is at least 3x the second largest
        if len(lengths) > 1 and lengths[0][1] >= 3 * lengths[1][1]:
            result = lengths[0][0]
        # Check if Higher Order Repeat
        if result == "other":
            hor_kmers = generateKmersFromFastaForward(seq, 6, True)
            hor_list = list(islice(hor_kmers, 1, size))
            x = top_n_frequent_distances(calculate_distances(hor_list),4)
            result = classify_hor(x)
            # Write to file
        new_annotations.append(f"{chrom}\t{chromStart}\t{chromEnd}\t{result}\t{score}\t{strand}\t{thickStart}\t{thickEnd}\t{get_color(result)}\n")
    with open(bed_file, "w") as bedfile:
        for entry2 in new_annotations:
        # Assuming entry is a tuple or list and you want to join the parts with tabs
            
            bedfile.write(entry2)

    print(f"Satellites annotated in {bed_file}!\n")


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
    return top_n

def get_unique_first_column_values(bed_file_path):
    unique_values = set()
    
    with open(bed_file_path, 'r') as file:
        next(file)  # Skip the first row
        for line in file:
            columns = line.strip().split('\t')
            if columns:
                unique_values.add(columns[0])
    
    return list(unique_values)

def classify_hor(distances_counts):
    if not isinstance(distances_counts, list) or not all(isinstance(i, tuple) and len(i) == 2 for i in distances_counts):
        raise ValueError("Input must be a list of (distance, count) tuples.")
    
    distances, counts = zip(*distances_counts)
    top_distance, top_count = distances[0], counts[0]
    
    if not (169 <= top_distance <= 173):
        #Check if multiple is larger
        if any(top_distance % x == 0 for x in (174, 173, 172, 171, 170, 169, 168)) and (top_distance > 174):
            return "active_hor"
        else:
            if top_distance > 4:
                return f"other_{top_distance}bp_repeat"
            else:
                return "other"
    
    half_top_count = top_count / 1.5
    multiples_within_half = False
    
    for i in range(1, len(distances)):
        if (169 <= distances[i] <= 173):
            continue
        modd = distances[i] % top_distance
        if modd in [0, top_distance-1, 1, top_distance-2, 2, top_distance-3, 3, top_distance-4, 4, top_distance-5, 5, top_distance-6, 6]:  # Allow small deviations
            if int(counts[i]) >= int(half_top_count):
                multiples_within_half = True
                break  # No need to check further if condition met
    
    return "active_hor" if multiples_within_half else "hor"

def get_color(label):
    return COLOR_MAPPING.get(label, COLOR_MAPPING['other']) if label in COLOR_MAPPING else (
        COLOR_MAPPING['other'] if label.startswith("other") else "0,0,0"
    )
