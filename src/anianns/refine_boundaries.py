import polars as pl
from anianns.kmer_utils import generate_kmers_from_fasta
from anianns.general_utils import extract_region
from itertools import islice
import math
import sys

from anianns.classification import load_all_kmer_dbs, classify_kmers

def detect_precise_boundaries(fasta_file, seq_id, seq_len, window, k, coordinates, verbose=False, classify=False, bordering=False):
    # Get the sequence of the specified region, subtracted by 1.5 * window size as a potential buffer region.
    interval = math.ceil(window / 2)
    array_seq_size = coordinates[1] - coordinates[0]
    
    # Check if the buffer is too large for the array size
    if array_seq_size < (2*window):
        print(f"Scaling down window size for {array_seq_size} bp array\n")
        window = math.ceil(array_seq_size // 4)
        interval = math.ceil(window // 2)

    core_start = coordinates[0] + window + interval
    core_end = coordinates[1] - window - interval
    core_seq_size = core_end - core_start

    if verbose:
        print(f"Estimated satellite array is {coordinates}, equal to {array_seq_size} bp\n")
    core_seq = extract_region(
        fasta_file=fasta_file,
        chr=seq_id,
        region_start = core_start,
        region_end= core_end
    )

    left_boundary = detect_left_boundary(
        fasta_file = fasta_file, 
        array_seq = core_seq,
        array_seq_size = core_seq_size,
        seq_id = seq_id, 
        boundary_point = coordinates[0],
        limit = None,
        k = k, 
        window = window,
        interval = interval,
        boundary_chunk_size=100,
        verbosity=verbose,
        bordering=bordering,
        boundary=0
    )

    right_boundary = detect_right_boundary(
        fasta_file = fasta_file, 
        array_seq = core_seq,
        array_seq_size = core_seq_size,
        seq_id = seq_id, 
        boundary_point = coordinates[1],
        limit = None,
        k = k, 
        window = window,
        interval = interval,
        boundary_chunk_size=100,
        verbosity=verbose,
        bordering=bordering,
        boundary=seq_len
    )

    if not core_seq:
        print(f"Unable to get boundary for {core_start}-{core_end}\n")
        return None
    
    if not left_boundary or not right_boundary:
        print(f"Error getting boundaries for {core_start}-{core_end}\n")
        return None

    if classify:
        array_seq_size = right_boundary - left_boundary
        array_seq = extract_region(
            fasta_file=fasta_file,
            chr=seq_id,
            region_start = left_boundary,
            region_end=right_boundary
        )
        if not array_seq:
            print(f"Unable to get boundary for {coordinates}\n")
            return None
        
        array_kmers = generate_kmers_from_fasta(array_seq, k, True)
        query_set = set(islice(array_kmers, array_seq_size))
        
        # Use the classification function
        if not verbose:
            verbose = False
        best_match, results = classify_kmers(query_set, classify, verbose)
        if not best_match:
            best_match = "Unknown"
        if verbose:
            print("---------------------")
        if verbose:
            if best_match:
                print(f"Best classification: {best_match}")
            else:
                print("No significant match found")

    

    return (left_boundary, right_boundary, best_match) if classify else (left_boundary, right_boundary, None)


def detect_left_boundary(fasta_file, array_seq, array_seq_size, seq_id, boundary_point, limit, k, window, interval, boundary_chunk_size, verbosity=False, bordering=False, boundary=0):
    
    boundary_start = boundary_point - window - interval
    # Ensure value is less than 0
    if boundary_start < boundary:
        boundary_start = boundary
    boundary_end = boundary_point + window + interval
    boundary_size = boundary_end - boundary_start
    border_seq = extract_region(
        fasta_file=fasta_file,
        chr=seq_id,
        region_start = boundary_start,
        region_end=boundary_end
    )
    # Generate k-mers and convert to sets for faster operations
    array_kmers = generate_kmers_from_fasta(array_seq, k, True)
    border_kmers = generate_kmers_from_fasta(border_seq, k, True)
    
    array_kmer_set = set(islice(array_kmers, array_seq_size))
    array_kmers = generate_kmers_from_fasta(array_seq, k, True)
    border_kmer_list = list(islice(border_kmers, boundary_size))

    #Process boundary in fixed length of 100bp
    steps = math.ceil(boundary_size / boundary_chunk_size)
    last_nonzero_step = find_target_fixed_window_left(
        array_kmer_set = array_kmer_set, 
        border_kmer_list = border_kmer_list, 
        boundary_size = boundary_size, 
        offset = boundary_start,
        verbose=verbosity,
        step_size = steps
    )
    possible_range = (boundary_start + (boundary_chunk_size * last_nonzero_step),boundary_start + (boundary_chunk_size * last_nonzero_step) + 2*boundary_chunk_size)

    # Case where boundary needs to be extended to the left
    if last_nonzero_step <= 1 :
        if verbosity:
            print("Boundary needs to be extended left. Likely an error\n")

        return None

    # Case where boundary needs to be extended to the right
    elif last_nonzero_step >= steps - 1:
        if verbosity:
            print("Boundary needs to be extended right.\n")
        updated_boundary_start = boundary_end - interval
        updated_boundary_end = boundary_end + (interval * 4)
        updated_boundary_size = updated_boundary_end - updated_boundary_start
        updated_border_seq = extract_region(
            fasta_file=fasta_file,
            chr=seq_id,
            region_start = updated_boundary_start,
            region_end=updated_boundary_end
        )
        updated_border_kmers = generate_kmers_from_fasta(updated_border_seq,k,True)
        updated_border_kmer_list = list(islice(updated_border_kmers,updated_boundary_size))
        updated_steps = math.ceil(updated_boundary_size / 100)
        last_nonzero_step = find_target_fixed_window_left(
            array_kmer_set = array_kmer_set, 
            border_kmer_list = updated_border_kmer_list, 
            boundary_size = updated_boundary_size, 
            offset = updated_boundary_start,
            verbose=verbosity,
            step_size = updated_steps
        )

        # Error case
        if last_nonzero_step <= 1 or last_nonzero_step >= updated_steps -1:
            if verbosity:
                print("Unable to resolve left boundary\n")
            return None
        
        #find_target_fixed_window_left(array_kmer_set, updated_border_kmer_list, updated_boundary_size, updated_boundary_start, updated_steps)
        possible_range = (updated_boundary_start + (100 * last_nonzero_step), updated_boundary_start + (100 * last_nonzero_step) + 200)
        #print(f"New possible range: {possible_range[0]}-{possible_range[1]}")
        if verbosity:
            print(f"Narrowing down the range for the updated left boundary: {possible_range[0]-math.ceil(boundary_chunk_size/2)}-{possible_range[1] + math.ceil(boundary_chunk_size/2)}\n")
        border_array_start = (boundary_chunk_size * last_nonzero_step)
        border_array_end = border_array_start + (2*boundary_chunk_size)

        estimated_index = find_last_matching_index(border_kmer_list, array_kmer_set, border_array_start, border_array_end)
        if verbosity:
            print(f"Estimated left boundary: {estimated_index + updated_boundary_start + k}\n") 
        return estimated_index + updated_boundary_start + k


    else:
        if verbosity:
            print(verbosity)
            print(f"Narrowing down the range for the left boundary: {possible_range[0]-math.ceil(boundary_chunk_size/2)}-{possible_range[1] + math.ceil(boundary_chunk_size/2)}\n")
        border_array_start = (boundary_chunk_size * last_nonzero_step)
        border_array_end = border_array_start + (2*boundary_chunk_size)

        estimated_index = find_last_matching_index(border_kmer_list, array_kmer_set, border_array_start, border_array_end)
        if verbosity:
            print(f"Estimated left boundary: {estimated_index + boundary_start + k}\n")
        return estimated_index + boundary_start - k  # Return the estimated boundary position adjusted by k-mer size

def detect_right_boundary(fasta_file, array_seq, array_seq_size, seq_id, boundary_point, limit, k, window, interval, boundary_chunk_size, verbosity=False, bordering=False, boundary=0):
    boundary_start = boundary_point - window - interval
    boundary_end = boundary_point + window + interval
    boundary_start = boundary_point - window - interval
    # Ensure value is less than 0
    if boundary_end > boundary:
        boundary_end = boundary
    boundary_size = boundary_end - boundary_start
    border_seq = extract_region(
        fasta_file=fasta_file,
        chr=seq_id,
        region_start = boundary_start,
        region_end=boundary_end
    )
    
    # Generate k-mers and convert to sets for faster operations
    array_kmers = generate_kmers_from_fasta(array_seq, k, True)
    border_kmers = generate_kmers_from_fasta(border_seq, k, True)
    
    array_kmer_set = set(islice(array_kmers, array_seq_size))
    border_kmer_list = list(islice(border_kmers, boundary_size))

    #Process boundary in fixed length of 100bp
    steps = math.ceil(boundary_size / boundary_chunk_size)
    last_nonzero_step = find_target_fixed_window_right(
        array_kmer_set = array_kmer_set, 
        border_kmer_list = border_kmer_list, 
        boundary_size = boundary_size, 
        offset = boundary_start,
        verbose=verbosity,
        step_size = steps
        )
    possible_range = (boundary_start + (boundary_chunk_size * last_nonzero_step),boundary_start + (boundary_chunk_size * last_nonzero_step) + 2*boundary_chunk_size)

    # Case for boundary extending to the right
    if last_nonzero_step >= steps -1:
        if verbosity:
            print("Boundary needs to be extended right. Likely an error.\n")
        return None

    # Case for boundary extending left
    elif last_nonzero_step <= 1:
        if verbosity:
            print("Boundary needs to be extended left.\n")
        updated_boundary_end = boundary_start + interval
        updated_boundary_start = boundary_start - (interval * 4)
        updated_boundary_size = updated_boundary_end - updated_boundary_start
        updated_border_seq = extract_region(
            fasta_file=fasta_file,
            chr=seq_id,
            region_start = updated_boundary_start,
            region_end=updated_boundary_end
        )
        updated_border_kmers = generate_kmers_from_fasta(updated_border_seq,k,True)
        updated_border_kmer_list = list(islice(updated_border_kmers,updated_boundary_size))
        updated_steps = math.ceil(updated_boundary_size / 100)
        last_nonzero_step = find_target_fixed_window_right(
            array_kmer_set = array_kmer_set, 
            border_kmer_list = updated_border_kmer_list, 
            boundary_size = updated_boundary_size, 
            offset = updated_boundary_start,
            verbose=verbosity,
            step_size = updated_steps
        )

        if verbosity:
            print(f"Narrowing down the range for the updated right boundary: {possible_range[0]-math.ceil(boundary_chunk_size/2)}-{possible_range[1] + math.ceil(boundary_chunk_size/2)}\n")
        border_array_start = (boundary_chunk_size * last_nonzero_step)
        border_array_end = border_array_start + (2*boundary_chunk_size)

        estimated_index = find_last_matching_index(border_kmer_list, array_kmer_set, border_array_start, border_array_end)
        if verbosity:
            print(updated_boundary_start)
            print(f"Estimated right boundary: {estimated_index + updated_boundary_start + k}\n") 
        return estimated_index + updated_boundary_start + k



    else:
        if verbosity:
            print(f"Narrowing down the range for the right boundary: {possible_range[0]-math.ceil(boundary_chunk_size/2)}-{possible_range[1] + math.ceil(boundary_chunk_size/2)}")
        border_array_start = (boundary_chunk_size * last_nonzero_step)
        border_array_end = border_array_start + (2*boundary_chunk_size)

        estimated_index = find_last_matching_index(border_kmer_list, array_kmer_set, border_array_start, border_array_end)
        if verbosity:
            print(border_array_start)
            print(f"Estimated right boundary: {estimated_index + border_array_start + k}\n")
        return estimated_index + boundary_start + k  # Return the estimated boundary position adjusted by k-mer size

def find_target_fixed_window_left(array_kmer_set, border_kmer_list, boundary_size, offset, verbose, step_size=100):
    steps = math.ceil(boundary_size / step_size)
    threshold = math.ceil(step_size / 4)
    last_nonzero_step = steps - 1

    shared_counts = []

    count = 0
    for i in range(step_size - 1, -1, -1):
        start_idx = i * steps
        end_idx = min((i + 1) * steps, len(border_kmer_list))
        border_chunk_set = set(border_kmer_list[start_idx:end_idx])
        shared_count = len(array_kmer_set & border_chunk_set)
        shared_counts.append((i, shared_count))

        if shared_count > threshold and count < 10:
            last_nonzero_step = i
            count = 0
        else:
            count += 1

    # Normalize for plotting
    max_count = max(count for _, count in shared_counts) or 1
    plot_height = 10  # Reduced from 20 to 10 for half the height
    scaled_counts = [
        (i, int((count / max_count) * plot_height)) for i, count in shared_counts
    ]

    if verbose:
        print(f"Find left boundary point from range {offset}-{offset + (steps * step_size)}\n")
        for y in range(plot_height, -1, -1):
            line = f"{str(int(y * max_count / plot_height)).rjust(4)} | "
            for idx, height in reversed(scaled_counts):
                if height >= y:
                    # Mark the stop position with a distinct symbol, e.g. '|'
                    if idx == last_nonzero_step:
                        line += "|"
                    else:
                        line += "#"
                else:
                    line += " "
            print(line)

        # X-axis line
        axis_line = "     " + "-" * step_size
        print(axis_line)

        # X-axis labels: leftmost and rightmost only
        left_label = str(offset)
        right_label = str(offset + (steps * step_size))
        spacing = step_size - len(left_label) - len(right_label)
        label_line = "     " + left_label + (" " * spacing) + right_label
        print(label_line)

        print(f"\nBoundary stop index: {last_nonzero_step}\n")

    # Return the index of the left boundary (slightly adjusted)
    return last_nonzero_step - 1

def find_target_fixed_window_right(array_kmer_set, border_kmer_list, boundary_size, offset, verbose, step_size=100):
    steps = math.ceil(boundary_size / step_size)
    threshold = math.ceil(step_size / 4)
    last_nonzero_step = 0

    shared_counts = []
    counter = 0
    for i in range(step_size):
        start_idx = i * steps
        end_idx = min((i + 1) * steps, len(border_kmer_list))
        border_chunk_set = set(border_kmer_list[start_idx:end_idx])
        shared_count = len(array_kmer_set & border_chunk_set)
        shared_counts.append((i, shared_count))
        
        if shared_count > threshold and counter < 10:
            last_nonzero_step = i
            counter = 0
        else:
            counter += 1

    # Normalize for plotting
    max_count = max(count for _, count in shared_counts) or 1
    plot_height = 10  # match left-hand version
    scaled_counts = [
        (i, int((count / max_count) * plot_height)) for i, count in shared_counts
    ]

    if verbose:
        print(f"Find right boundary point from range {offset}-{offset + (steps * step_size)}\n")

        for y in range(plot_height, -1, -1):
            line = f"{str(int(y * max_count / plot_height)).rjust(4)} | "
            for idx, height in scaled_counts:
                if height >= y:
                    # Highlight the boundary stop position with '|'
                    if idx == last_nonzero_step:
                        line += "|"
                    else:
                        line += "#"
                else:
                    line += " "
            print(line)

        # X-axis line
        axis_line = "     " + "-" * step_size
        print(axis_line)

        # X-axis labels: leftmost and rightmost
        left_label = str(offset)
        right_label = str(offset + (steps * step_size))
        spacing = step_size - len(left_label) - len(right_label)
        label_line = "     " + left_label + (" " * spacing) + right_label
        print(label_line)

        print(f"\nBoundary stop index: {last_nonzero_step}\n")

    return last_nonzero_step


def find_last_matching_index(border_kmer_list, array_kmer_set, border_array_start, border_array_end):
    try:
        for i in range(border_array_end - 1, border_array_start - 1, -1):
            kmer = border_kmer_list[i]
            if kmer in array_kmer_set:
                return i
    except IndexError:
        # TODO FIX!!
        #print("IndexError: border_array_start or border_array_end is out of bounds.")
        # If we hit an index error, we can assume the border array start is the best estimate
        # This is a fallback to ensure we return a valid index
        #print(f"Returning border_array_start: {border_array_start},{border_array_end}")
        pass
    return border_array_start



def report_borders(fa, seq_id, seq_len, band, offset, window, k, df: pl.DataFrame, classify, verbose, quiet) -> None:
    band = int(band) * 1_000_000
    if not quiet:
        print(f"Inferring satellite locations & boundaries for {seq_id}:\n")
    # Load k-mer sets here for classification
    loaded_kmer_dbs = None
    if classify:

        loaded_kmer_dbs = load_all_kmer_dbs(classify)
        if verbose:
            print(f"Loaded {len(loaded_kmer_dbs)} k-mer database(s) for classification.")
        # Or iterate through all databases and create supersets
        supersets_dict = {}
        for db_name, (k_val, sets_dict) in loaded_kmer_dbs.items():
            #print(f"Database {db_name}: k={k_val}, sets={list(sets_dict.keys())}")
            
            # Create a super_set containing all k-mers from all sets in this database
            super_set = set()
            for set_name, kmers in sets_dict.items():
                super_set.update(kmers)
            
            # Store the superset in the dictionary
            supersets_dict[db_name] = super_set
            if verbose:
                print(f"Loaded {db_name}: {len(super_set)} total k-mers")

    else:
        supersets_dict = False
        
        # Now supersets_dict contains: {"db_name": {combined_kmers}, ...}
    
    # pull out start/end as plain Python lists for easy indexing
    if verbose:
        print("\n")
        print(f"Processing {seq_id} boundaries with window size {window} and k-mer size {k}:\n")
    starts = [s - offset for s in df["start"].to_list()]
    ends   = [e - offset for e in df["end"].to_list()]
    # Ensure all values are greater than 0
    starts = [s if s > 0 else 0 for s in starts]
    ends = [e if e > 0 else 0 for e in ends]
    new_starts = []
    new_ends = []
    classification = []
    n = df.height
    assert len(starts) == len(ends)
    i = 0
    while i < len(starts):
        next_i = i + 1
        potential = ends[i] % band

        if potential <= window*3 or potential >= band - (window*3):
            if verbose:
                print(f"Band border region detected at end position {ends[i]}")
            if next_i < len(starts):
                ends[i] = ends[next_i]
                del starts[next_i]
                del ends[next_i]
            coordinates = (int(starts[i]), int(ends[i]))
            try:
                updated_boundaries = detect_precise_boundaries(fa, seq_id, seq_len, window, k, coordinates, verbose, supersets_dict)
                if verbose:
                    print(f"Estimated new boundaries: {updated_boundaries[0]}-{updated_boundaries[1]}")
                new_starts.append(updated_boundaries[0])
                new_ends.append(updated_boundaries[1])
                classification.append(updated_boundaries[2])
            except Exception as e:
                print(f"Error occurred while updating boundaries: {e}")
            # Do not increment i, as the next region is now at the same index

        else:

            # check if end of current row == start of next
            if next_i < len(starts) and abs(ends[i] - starts[next_i]) <= 4 * window:
                coordinates = (int(starts[i]), int(ends[i]))
                try:
                    updated_boundaries = detect_precise_boundaries(fa, seq_id, seq_len, window, k, coordinates, verbose, supersets_dict)
                    if verbose:
                        print(f"Estimated new boundaries: {updated_boundaries[0]}-{updated_boundaries[1]}")
                    new_starts.append(updated_boundaries[0])
                    new_ends.append(updated_boundaries[1])
                    classification.append(updated_boundaries[2])
                except Exception as e:
                    print(f"Error occurred while updating boundaries: {e}")
            else:
                coordinates = (int(starts[i]), int(ends[i]))
                updated_boundaries = detect_precise_boundaries(fa, seq_id, seq_len, window, k, coordinates, verbose, supersets_dict)
                if not updated_boundaries:
                    new_starts.append(0)
                    new_ends.append(1)
                    classification.append("Unknown")
                else:
                    new_starts.append(updated_boundaries[0])
                    new_ends.append(updated_boundaries[1])
                    classification.append(updated_boundaries[2])
            if verbose:
                print("--------------------------------------------------\n")
        i += 1
    return new_starts, new_ends, classification