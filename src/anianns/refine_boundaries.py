import polars as pl
from anianns.kmer_utils import (
    generate_kmers_from_fasta,
    generate_kmers_from_fasta_forward_only,
)
from itertools import islice
import math
import sys
import numpy as np

from anianns.classification import load_all_kmer_dbs, classify_kmers

from anianns.general_utils import (
    extract_region,
    extract_histograms_by_name,
    generate_kmers_from_fasta,
    generate_kmers_from_fasta_forward_only,
    top_n_frequent_distances,
    calculate_distances,
    merge_close_values,
)


def refine_L_from_peaks(dist_counts, L_init, max_rel_err=0.1, max_n=40):
    """
    Refine monomer length L from observed distances using a simple
    regression: distance_i ≈ n_i * L.

    dist_counts : list of (distance, value)
    L_init      : initial guess of L (e.g. 171)
    max_rel_err : max relative error |d - n*L_init| / (n*L_init) to
                  accept a peak as a harmonic when building the model
    max_n       : maximum harmonic index to consider
    """
    distances = np.array([d for d, v in dist_counts])

    ns = []
    ds = []

    for d in distances:
        if d <= 0:
            continue
        n = int(round(d / L_init))
        if n < 1 or n > max_n:
            continue
        rel_err = abs(d - n * L_init) / (n * L_init)
        if rel_err <= max_rel_err:
            ns.append(n)
            ds.append(d)

    if len(ns) < 2:
        # not enough data to refine; fall back
        return L_init

    ns = np.array(ns, dtype=float)
    ds = np.array(ds, dtype=float)

    # Fit slope in d ≈ n * L (through origin)
    L_refined = np.sum(ns * ds) / np.sum(ns**2)
    return L_refined


def collect_harmonics_with_rel_tol(
    dist_counts, L, max_harmonics=40, rel_tol=0.02, abs_tol_min=5.0
):
    """
    For each harmonic index n, find the strongest peak near n*L,
    allowing a relative tolerance (and a minimum absolute tolerance).

    rel_tol     : allowed relative deviation |d - nL| / (nL)
    abs_tol_min : minimum absolute tolerance in bp
    """
    dist_counts = sorted(dist_counts, key=lambda x: x[0])
    distances = np.array([d for d, v in dist_counts])
    values = np.array([v for d, v in dist_counts])

    harmonic_ns = []
    harmonic_vals = []
    harmonic_ds = []

    for n in range(1, max_harmonics + 1):
        target = n * L
        abs_tol = max(abs_tol_min, rel_tol * target)

        mask = (distances >= target - abs_tol) & (distances <= target + abs_tol)
        if np.any(mask):
            idx_local = np.argmax(values[mask])
            val = values[mask][idx_local]
            d = distances[mask][idx_local]

            harmonic_ns.append(n)
            harmonic_vals.append(val)
            harmonic_ds.append(d)

    return (
        np.array(harmonic_ns, dtype=int),
        np.array(harmonic_ds, dtype=float),
        np.array(harmonic_vals, dtype=float),
    )


def hor_test_local_enrichment(
    dist_counts,
    L_init,
    refine=True,
    max_rel_err_refine=0.1,
    max_harmonics=40,
    rel_tol=0.02,
    abs_tol_min=5.0,
    enrichment_threshold=2.0,
    min_n_for_HOR=2,
):
    """
    Detect HORs using locally enriched harmonics with approximate spacing.

    dist_counts : list of (distance, value) (counts or proportions)
    L_init      : initial guess of monomer length
    refine      : whether to refine L using smaller harmonics
    enrichment_threshold : how much bigger than neighbors a harmonic
                           must be to be considered HOR (e.g. 2.0 = 2x)
    min_n_for_HOR        : only consider harmonics with n >= this as HOR
    """
    # 1. refine L (optional)
    if refine:
        L = refine_L_from_peaks(
            dist_counts, L_init, max_rel_err=max_rel_err_refine, max_n=max_harmonics
        )
    else:
        L = float(L_init)

    # 2. collect harmonics with relative tolerance
    ns, ds, H = collect_harmonics_with_rel_tol(
        dist_counts,
        L,
        max_harmonics=max_harmonics,
        rel_tol=rel_tol,
        abs_tol_min=abs_tol_min,
    )

    enrichments = []
    hor_flags = []

    for i, n in enumerate(ns):
        left = H[i - 1] if i - 1 >= 0 else None
        right = H[i + 1] if i + 1 < len(H) else None
        neighbors = [x for x in (left, right) if x is not None]

        if not neighbors:
            enrichments.append(np.nan)
            hor_flags.append(False)
            continue

        baseline = np.mean(neighbors)
        if baseline <= 0:
            enrichments.append(np.nan)
            hor_flags.append(False)
            continue

        ratio = H[i] / baseline
        enrichments.append(ratio)

        hor_flags.append((n >= min_n_for_HOR) and (ratio >= enrichment_threshold))

    enrichments = np.array(enrichments, dtype=float)
    hor_flags = np.array(hor_flags, dtype=bool)

    hor_harmonics = list(ns[hor_flags])

    if len(hor_harmonics) == 0:
        return False
    else:
        return True, hor_harmonics


def ntr_prism(region, size, kmer_large, verbose=False):
    kmers_low = generate_kmers_from_fasta_forward_only(seq=region, k=6, quiet=True)
    kmer_list_low = list(islice(kmers_low, 1, size))
    histo_list_low = calculate_distances(kmer_list_low)
    top_dists_low = top_n_frequent_distances(histo_list_low, 6)
    grouped_dists_low = merge_close_values(top_dists_low, 6)

    # Defaults so we never return unbound variables
    monomer_size = None
    is_hor = False
    periodicity = None
    nested_repeats = []

    if size <= 0:
        if verbose:
            print("Size must be > 0 to compute percentages.")
        return None, False, None

    grouped = list(grouped_dists_low)
    if not grouped:
        if verbose:
            print("No grouped distances found.")
        return None, False, None

    top_total_count = sum(count for (_, count) in grouped)

    if verbose:
        print("Top distance groups (distance: count, percentage of size):")
    for dist, count in grouped:
        pct = (count / size) * 100
        if verbose:
            print(f"  {dist}: {count} ({pct:.2f}%)")

    overall_pct = (top_total_count / size) * 100
    if verbose:
        print(
            f"\nTotal in top {len(grouped)} distances: {top_total_count} / {size} = {overall_pct:.2f}%\n"
        )

    if overall_pct < 10.0:
        if verbose:
            print(
                "Warning: Top distance groups account for less than 10% of the total size. Removing satellite.\n"
            )
        return None, False, None

    if verbose:
        print(f"low: {grouped}\n")

    top10 = grouped[:10]
    if not top10:
        return None, False, None

    # Lowest distance by value
    lowest_dist, lowest_count = min(top10, key=lambda x: x[0])

    # Most populous by count
    most_populous_dist, most_populous_count = max(top10, key=lambda x: x[1])

    # Decide monomer_size
    monomer_size = lowest_dist  # sensible default

    if lowest_dist is None or lowest_dist == 0:
        if verbose:
            print("Lowest distance is 0/None; cannot check harmonics safely.")
        monomer_size = most_populous_dist  # fallback to most populous
    elif lowest_dist != most_populous_dist:
        ratio = most_populous_dist / lowest_dist
        n = round(ratio)

        rel_tol = 0.02  # 2% tolerance

        # NOTE: compare absolute error in ratio space
        is_harmonic = (n >= 1) and (abs(ratio - n) <= rel_tol)

        if is_harmonic:
            if verbose:
                print(f"harmonic: {most_populous_dist} ≈ {n} × {lowest_dist}")
            monomer_size = lowest_dist
        else:
            if verbose:
                print("Not a harmonic relationship")
            # keep monomer_size = lowest_dist (already set)

    # HOR test — uses monomer_size as L_init (not always lowest_dist)
    res = hor_test_local_enrichment(
        grouped,
        L_init=monomer_size,
        refine=True,
        max_harmonics=400,
        rel_tol=0.02,
        abs_tol_min=5.0,
        enrichment_threshold=3.0,
        min_n_for_HOR=2,
    )

    # hor_test_local_enrichment returns False OR (True, hor_harmonics)
    if not res:
        if verbose:
            print("No HOR detected based on local enrichment of harmonics.\n")
    else:
        # res is (True, hor_harmonics)
        periodicity = res[1]
        # FInd the value that matches the hor_harmonics
        if verbose:
            print(periodicity)

        is_hor = True
        tolerance = 0.02  # 2%
        for i in periodicity:
            target = monomer_size * i
            match = next(
                (t for t in top10 if abs(t[0] - target) <= 0.02 * target), None
            )

            nested_repeats.append(match[0])

    return monomer_size, is_hor, periodicity, nested_repeats


def detect_precise_boundaries(
    fasta_file,
    seq_id,
    seq_len,
    window,
    k,
    coordinates,
    verbose,
    classify,
    previous_coordinates,
    interval=None,
):
    # Get the sequence of the specified region, subtracted by 1.5 * window size as a potential buffer region.
    interval = math.ceil(window / 2)
    array_seq_size = coordinates[1] - coordinates[0]

    # Check if the buffer is too large for the array size
    if array_seq_size < (2 * window):
        print(f"Scaling down window size for {array_seq_size} bp array\n")
        window = math.ceil(array_seq_size // 4)
        interval = math.ceil(window // 2)

    core_start = coordinates[0] + window + interval
    core_end = coordinates[1] - window - interval
    core_seq_size = core_end - core_start

    if verbose:
        print(
            f"Estimated satellite array is {coordinates}, equal to {array_seq_size} bp\n"
        )

    core_seq = extract_region(
        fasta_file=fasta_file, chr=seq_id, region_start=core_start, region_end=core_end
    )

    if not core_seq:
        if verbose:
            print(f"Unable to get core sequence for {core_start}-{core_end}\n")
        return None

    prism_res = ntr_prism(core_seq, core_seq_size, k)

    # If ntr_prism says "remove satellite", drop this whole region
    if prism_res[0] is None:
        if verbose:
            print(f"ntr_prism returned None for {coordinates}; removing region.\n")
        return None

    # This survived the ntr_prism purge, find the boundary
    left_boundary = detect_left_boundary(
        fasta_file=fasta_file,
        array_seq=core_seq,
        array_seq_size=core_seq_size,
        seq_id=seq_id,
        boundary_point=coordinates[0],
        k=k,
        window=window,
        interval=interval,
        boundary_chunk_size=100,
        verbosity=verbose,
        prev_border_coordinate=previous_coordinates[1],
        boundary=0,
        expanded=False,
    )

    right_boundary = detect_right_boundary(
        fasta_file=fasta_file,
        array_seq=core_seq,
        array_seq_size=core_seq_size,
        seq_id=seq_id,
        boundary_point=coordinates[1],
        limit=None,
        k=k,
        window=window,
        interval=interval,
        boundary_chunk_size=100,
        verbosity=verbose,
        bordering=False,
        boundary=seq_len,
    )

    if not left_boundary or not right_boundary:
        if verbose:
            print(f"Error getting boundaries for {core_start}-{core_end}\n")
        return None

    if classify:
        array_seq_size = right_boundary - left_boundary
        array_seq = extract_region(
            fasta_file=fasta_file,
            chr=seq_id,
            region_start=left_boundary,
            region_end=right_boundary,
        )
        if not array_seq:
            if verbose:
                print(f"Unable to extract array sequence for {coordinates}\n")
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

    return (
        (left_boundary, right_boundary, best_match, prism_res)
        if classify
        else (left_boundary, right_boundary, None, prism_res)
    )


def detect_left_boundary(
    fasta_file,
    array_seq,
    array_seq_size,
    seq_id,
    boundary_point,
    k,
    window,
    interval,
    boundary_chunk_size,
    verbosity,
    prev_border_coordinate,
    boundary,
    expanded=False,
):
    boundary_start = boundary_point - window - (interval * 4)
    # if prev_border_coordinate > (boundary_start - interval):
    # print(f"Potential conflict as previous boundary at {prev_border_coordinate} close to {boundary_start}")

    # Ensure the boundary start isn't less than 0
    if boundary_start < boundary:
        boundary_start = boundary

    boundary_end = boundary_point + window + (interval * 4)
    boundary_size = boundary_end - boundary_start

    border_seq = extract_region(
        fasta_file=fasta_file,
        chr=seq_id,
        region_start=boundary_start,
        region_end=boundary_end,
    )

    # Generate k-mers and convert to sets for faster operations
    array_kmers = generate_kmers_from_fasta(array_seq, k, True)
    border_kmers = generate_kmers_from_fasta(border_seq, k, True)

    array_kmer_set = set(islice(array_kmers, array_seq_size))
    array_kmers = generate_kmers_from_fasta(array_seq, k, True)
    border_kmer_list = list(islice(border_kmers, boundary_size))

    # Process boundary in fixed length of 100bp
    steps = math.ceil(boundary_size / boundary_chunk_size)
    last_nonzero_step = find_target_fixed_window_left(
        array_kmer_set=array_kmer_set,
        border_kmer_list=border_kmer_list,
        boundary_size=boundary_size,
        offset=boundary_start,
        verbose=verbosity,
        step_size=steps,
    )
    possible_range = (
        boundary_start + (boundary_chunk_size * last_nonzero_step),
        boundary_start
        + (boundary_chunk_size * last_nonzero_step)
        + 2 * boundary_chunk_size,
    )

    # Case where boundary needs to be extended to the left
    if last_nonzero_step <= 1:
        if boundary_start <= window:
            return boundary_start
        else:
            if verbosity:
                print("Boundary needs to be extended left...\n")

            if not expanded:
                return detect_left_boundary(
                    fasta_file=fasta_file,
                    array_seq=array_seq,
                    array_seq_size=array_seq_size,
                    seq_id=seq_id,
                    boundary_point=boundary_start,  # This gets modified
                    k=k,
                    window=window,
                    interval=interval,
                    boundary_chunk_size=100,
                    verbosity=verbosity,
                    prev_border_coordinate=prev_border_coordinate,
                    boundary=0,
                    expanded=True,
                )
            else:
                if prev_border_coordinate > boundary_point - window:
                    prev_border_coordinate + 1
                    return prev_border_coordinate + 1

            """updated_boundary_start = int(boundary_start - (interval * 10))
            updated_boundary_end = int(boundary_end + (interval * 10))
            print(f"DDD {updated_boundary_start} - {updated_boundary_end}\n")
            updated_boundary_size = updated_boundary_end - updated_boundary_start
            updated_border_seq = extract_region(
                fasta_file=fasta_file,
                chr=seq_id,
                region_start=updated_boundary_start,
                region_end=updated_boundary_end,
            )
            updated_border_kmers = generate_kmers_from_fasta(updated_border_seq, k, True)
            updated_border_kmer_list = list(
                islice(updated_border_kmers, updated_boundary_size)
            )
            updated_steps = math.ceil(updated_boundary_size / 100)
            last_nonzero_step = find_target_fixed_window_left(
                array_kmer_set=array_kmer_set,
                border_kmer_list=updated_border_kmer_list,
                boundary_size=updated_boundary_size,
                offset=updated_boundary_start,
                verbose=verbosity,
                step_size=updated_steps,
            )
            print("Testing") """

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
            region_start=updated_boundary_start,
            region_end=updated_boundary_end,
        )
        updated_border_kmers = generate_kmers_from_fasta(updated_border_seq, k, True)
        updated_border_kmer_list = list(
            islice(updated_border_kmers, updated_boundary_size)
        )
        updated_steps = math.ceil(updated_boundary_size / 100)
        last_nonzero_step = find_target_fixed_window_left(
            array_kmer_set=array_kmer_set,
            border_kmer_list=updated_border_kmer_list,
            boundary_size=updated_boundary_size,
            offset=updated_boundary_start,
            verbose=verbosity,
            step_size=updated_steps,
        )

        # Error case
        if last_nonzero_step <= 1 or last_nonzero_step >= updated_steps - 1:
            if verbosity:
                print("Unable to resolve left boundary\n")
            return None

        # find_target_fixed_window_left(array_kmer_set, updated_border_kmer_list, updated_boundary_size, updated_boundary_start, updated_steps)
        possible_range = (
            updated_boundary_start + (100 * last_nonzero_step),
            updated_boundary_start + (100 * last_nonzero_step) + 200,
        )
        # print(f"New possible range: {possible_range[0]}-{possible_range[1]}")
        if verbosity:
            print(
                f"Narrowing down the range for the updated left boundary: {possible_range[0]-math.ceil(boundary_chunk_size/2)}-{possible_range[1] + math.ceil(boundary_chunk_size/2)}\n"
            )
        border_array_start = boundary_chunk_size * last_nonzero_step
        border_array_end = border_array_start + (2 * boundary_chunk_size)

        estimated_index = find_last_matching_index(
            border_kmer_list, array_kmer_set, border_array_start, border_array_end
        )
        if verbosity:
            print(
                f"Estimated left boundary: {estimated_index + updated_boundary_start + k}\n"
            )
        return estimated_index + updated_boundary_start + k

    else:
        if verbosity:
            print(verbosity)
            print(
                f"Narrowing down the range for the left boundary: {possible_range[0]-math.ceil(boundary_chunk_size/2)}-{possible_range[1] + math.ceil(boundary_chunk_size/2)}\n"
            )
        border_array_start = boundary_chunk_size * last_nonzero_step
        border_array_end = border_array_start + (2 * boundary_chunk_size)

        estimated_index = find_last_matching_index(
            border_kmer_list, array_kmer_set, border_array_start, border_array_end
        )
        if verbosity:
            print(f"Estimated left boundary: {estimated_index + boundary_start + k}\n")
        return (
            estimated_index + boundary_start - k
        )  # Return the estimated boundary position adjusted by k-mer size


def detect_right_boundary(
    fasta_file,
    array_seq,
    array_seq_size,
    seq_id,
    boundary_point,
    limit,
    k,
    window,
    interval,
    boundary_chunk_size,
    verbosity=False,
    bordering=False,
    boundary=0,
):
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
        region_start=boundary_start,
        region_end=boundary_end,
    )

    # Generate k-mers and convert to sets for faster operations
    array_kmers = generate_kmers_from_fasta(array_seq, k, True)
    border_kmers = generate_kmers_from_fasta(border_seq, k, True)

    array_kmer_set = set(islice(array_kmers, array_seq_size))
    border_kmer_list = list(islice(border_kmers, boundary_size))

    # Process boundary in fixed length of 100bp
    steps = math.ceil(boundary_size / boundary_chunk_size)
    last_nonzero_step = find_target_fixed_window_right(
        array_kmer_set=array_kmer_set,
        border_kmer_list=border_kmer_list,
        boundary_size=boundary_size,
        offset=boundary_start,
        verbose=verbosity,
        step_size=steps,
    )
    possible_range = (
        boundary_start + (boundary_chunk_size * last_nonzero_step),
        boundary_start
        + (boundary_chunk_size * last_nonzero_step)
        + 2 * boundary_chunk_size,
    )

    # Case for boundary extending to the right
    if last_nonzero_step >= steps - 1:
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
            region_start=updated_boundary_start,
            region_end=updated_boundary_end,
        )
        updated_border_kmers = generate_kmers_from_fasta(updated_border_seq, k, True)
        updated_border_kmer_list = list(
            islice(updated_border_kmers, updated_boundary_size)
        )
        updated_steps = math.ceil(updated_boundary_size / 100)
        last_nonzero_step = find_target_fixed_window_right(
            array_kmer_set=array_kmer_set,
            border_kmer_list=updated_border_kmer_list,
            boundary_size=updated_boundary_size,
            offset=updated_boundary_start,
            verbose=verbosity,
            step_size=updated_steps,
        )

        if verbosity:
            print(
                f"Narrowing down the range for the updated right boundary: {possible_range[0]-math.ceil(boundary_chunk_size/2)}-{possible_range[1] + math.ceil(boundary_chunk_size/2)}\n"
            )
        border_array_start = boundary_chunk_size * last_nonzero_step
        border_array_end = border_array_start + (2 * boundary_chunk_size)

        estimated_index = find_last_matching_index(
            border_kmer_list, array_kmer_set, border_array_start, border_array_end
        )
        if verbosity:
            print(updated_boundary_start)
            print(
                f"Estimated right boundary: {estimated_index + updated_boundary_start + k}\n"
            )
        return estimated_index + updated_boundary_start + k

    else:
        if verbosity:
            print(
                f"Narrowing down the range for the right boundary: {possible_range[0]-math.ceil(boundary_chunk_size/2)}-{possible_range[1] + math.ceil(boundary_chunk_size/2)}"
            )
        border_array_start = boundary_chunk_size * last_nonzero_step
        border_array_end = border_array_start + (2 * boundary_chunk_size)

        estimated_index = find_last_matching_index(
            border_kmer_list, array_kmer_set, border_array_start, border_array_end
        )
        if verbosity:
            print(border_array_start)
            print(
                f"Estimated right boundary: {estimated_index + border_array_start + k}\n"
            )
        return (
            estimated_index + boundary_start + k
        )  # Return the estimated boundary position adjusted by k-mer size


def find_target_fixed_window_left(
    array_kmer_set, border_kmer_list, boundary_size, offset, verbose, step_size=100
):
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
        print(
            f"Find left boundary point from range {offset}-{offset + (steps * step_size)}\n"
        )
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


def find_target_fixed_window_right(
    array_kmer_set, border_kmer_list, boundary_size, offset, verbose, step_size=100
):
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
        print(
            f"Find right boundary point from range {offset}-{offset + (steps * step_size)}\n"
        )

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

        print(
            f"\nBoundary stop index: {last_nonzero_step}, {offset + (step_size * last_nonzero_step)}\n"
        )

    return last_nonzero_step


def find_last_matching_index(
    border_kmer_list, array_kmer_set, border_array_start, border_array_end
):
    try:
        for i in range(border_array_end - 1, border_array_start - 1, -1):
            kmer = border_kmer_list[i]
            if kmer in array_kmer_set:
                return i
    except IndexError:
        # TODO FIX!!
        # print("IndexError: border_array_start or border_array_end is out of bounds.")
        # If we hit an index error, we can assume the border array start is the best estimate
        # This is a fallback to ensure we return a valid index
        # print(f"Returning border_array_start: {border_array_start},{border_array_end}")
        pass
    return border_array_start


def report_borders(
    fa: str,
    seq_id: str,
    seq_len: int,
    band: int,
    offset: int,
    window: int,
    k: int,
    df: pl.DataFrame,
    classify: str | bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """
    Infer satellite locations & precise boundaries for regions described in `df`.
    """
    # Convert band (in millions) to base-pairs
    band = int(band) * 1_000_000

    if not quiet:
        print(f"Inferring satellite locations & boundaries for {seq_id}...\n")

    # Load k-mer DBs (for classification) and build supersets (supersets_dict = classification dictionary)
    loaded_kmer_dbs = None
    if classify:
        loaded_kmer_dbs = load_all_kmer_dbs(classify)

        supersets_dict = {}
        for db_name, (k_val, sets_dict) in loaded_kmer_dbs.items():
            super_set = set()
            for set_name, kmers in sets_dict.items():
                super_set.update(kmers)
            supersets_dict[db_name] = super_set

            if verbose:
                print(f"Loaded {db_name}: {len(super_set)} total k-mers\n")
    else:
        supersets_dict = False
        if not quiet:
            print("No k-mer database provided.\n")

    if verbose:
        print("\n")
        print(
            f"Processing {seq_id} boundaries with window size {window} and k-mer size {k}:\n"
        )

    # Prepare starts/ends relative to the provided offset
    starts = [s - offset for s in df["start"].to_list()]
    ends = [e - offset for e in df["end"].to_list()]

    # Ensure non-negative
    starts = [s if s > 0 else 0 for s in starts]
    ends = [e if e > 0 else 0 for e in ends]

    new_starts = []
    new_ends = []
    classification = []
    monomer = []
    periodicity = []
    hor = []
    previous_new_boundaries = []

    assert len(starts) == len(ends)

    i = 0
    while i < len(starts):
        next_i = i + 1
        potential = ends[i] % band
        previous_coordinates = (
            previous_new_boundaries[-1] if previous_new_boundaries else (0, 1)
        )

        # Band-border / wrap-around handling
        # This code runs when the satellite boundary is by a band border
        if potential <= window * 3 or potential >= band - (window * 3):
            if verbose:
                print(f"Estimated borders cross a border region: {ends[i]+1}\n")
            if next_i < len(starts):
                # merge with next region
                ends[i] = ends[next_i]
                del starts[next_i]
                del ends[next_i]

            coordinates = (int(starts[i]) + 1, int(ends[i]) + 1)
            try:
                updated_boundaries = detect_precise_boundaries(
                    fasta_file=fa,
                    seq_id=seq_id,
                    seq_len=seq_len,
                    window=window,
                    k=k,
                    coordinates=coordinates,
                    verbose=verbose,
                    classify=supersets_dict,
                    previous_coordinates=previous_coordinates,
                )
                if verbose:
                    print(
                        f"Estimated new boundaries: {updated_boundaries[0]}-{updated_boundaries[1]}"
                    )
                # Various appendages
                new_starts.append(updated_boundaries[0])
                new_ends.append(updated_boundaries[1])
                classification.append(updated_boundaries[2])
                monomer.append(updated_boundaries[3][0])
                periodicity.append(updated_boundaries[3][2])
                hor.append(updated_boundaries[3][3])
                previous_new_boundaries.append(
                    (updated_boundaries[0], updated_boundaries[1])
                )
            except Exception as e:
                if verbose:
                    print(
                        f"Error occurred while updating boundaries at {coordinates} (merged border region): {e}"
                    )
            # Do not increment i, as the next region is now at the same index

        # This code runs if not by a band border
        else:
            coordinates = (int(starts[i]) + 1, int(ends[i]) + 1)
            try:
                updated_boundaries = detect_precise_boundaries(
                    fasta_file=fa,
                    seq_id=seq_id,
                    seq_len=seq_len,
                    window=window,
                    k=k,
                    coordinates=coordinates,
                    verbose=verbose,
                    classify=supersets_dict,
                    previous_coordinates=previous_coordinates,
                )
                if verbose:
                    print(
                        f"Estimated new boundaries: {updated_boundaries[0]}-{updated_boundaries[1]}"
                    )
                # Various appendages
                new_starts.append(updated_boundaries[0])
                new_ends.append(updated_boundaries[1])
                classification.append(updated_boundaries[2])
                monomer.append(updated_boundaries[3][0])
                periodicity.append(updated_boundaries[3][2])
                hor.append(updated_boundaries[3][3])
                previous_new_boundaries.append(
                    (updated_boundaries[0], updated_boundaries[1])
                )
            except Exception as e:
                if verbose:
                    print(
                        f"Error occurred while updating boundaries at {coordinates} (non-border region): {e}"
                    )

            if verbose:
                print("--------------------------------------------------\n")

        i += 1

    return new_starts, new_ends, classification, monomer, periodicity, hor
