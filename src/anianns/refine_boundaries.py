import polars as pl
import math
import numpy as np

from anianns.classification import load_all_kmer_dbs, classify_kmers

from anianns.general_utils import (
    extract_region,
    top_n_frequent_distances,
    merge_close_values,
)
from anianns.kmer_pipeline import canonical_kmer_hashes, forward_kmer_hashes
from anianns.kmer_utils import calculate_hash_distances


BOUNDARY_GRAPH_MAX_WIDTH = 80
BOUNDARY_GRAPH_PREFIX_WIDTH = 7


def boundary_graph_columns(
    shared_counts,
    stop_index,
    plot_height,
    max_columns=BOUNDARY_GRAPH_MAX_WIDTH - BOUNDARY_GRAPH_PREFIX_WIDTH,
):
    """Aggregate boundary-search chunks into terminal-safe graph columns."""
    if not shared_counts:
        return [], 1
    column_count = min(len(shared_counts), max(1, int(max_columns)))
    max_count = max(count for _, count in shared_counts) or 1
    columns = []
    for column in range(column_count):
        start = (column * len(shared_counts)) // column_count
        end = ((column + 1) * len(shared_counts)) // column_count
        bucket = shared_counts[start:end]
        bucket_count = max(count for _, count in bucket)
        columns.append(
            (
                any(index == stop_index for index, _ in bucket),
                int((bucket_count / max_count) * plot_height),
            )
        )
    return columns, max_count


def boundary_graph_label_line(offset, boundary_size, graph_width):
    """Format the coordinate labels without exceeding the graph width."""
    available_width = max(1, BOUNDARY_GRAPH_MAX_WIDTH - 5)
    graph_width = min(int(graph_width), available_width)
    left_label = str(offset)
    right_label = str(offset + boundary_size)
    if len(left_label) + len(right_label) >= graph_width:
        label_budget = max(1, (graph_width - 1) // 2)
        left_label = left_label[:label_budget]
        right_label = right_label[-label_budget:]
    spacing = max(1, graph_width - len(left_label) - len(right_label))
    return ("     " + left_label + (" " * spacing) + right_label)[
        :BOUNDARY_GRAPH_MAX_WIDTH
    ]


def recurrent_kmer_set(kmer_hashes, minimum_occurrences=3, minimum_unique=25):
    """Return repeat-specific k-mers suitable for a failed boundary search."""
    hashes = np.asarray(kmer_hashes, dtype=np.int32)
    hashes = hashes[hashes != 0]
    if hashes.size == 0:
        return None
    values, counts = np.unique(hashes, return_counts=True)
    recurrent = values[counts >= minimum_occurrences]
    if len(recurrent) < minimum_unique:
        return None
    return set(recurrent.tolist())


def format_ntrprism_verbose_value(kmer, prism_result, *, user_selected=False):
    """Format one NTRPrism result for an accepted-candidate summary."""
    label = f"NTRPrism k={int(kmer)}"
    if user_selected:
        label += " (user-selected)"
    monomer = prism_result[0]
    if monomer in (None, 0):
        return f"{label}: rejected"
    summary = f"{label}: monomer={monomer} bp"
    if prism_result[1]:
        summary += f"; HOR periodicity={prism_result[2]}"
    return summary


def region_canonical_hashes(
    fasta_file,
    seq_id,
    region_start,
    region_end,
    k,
    sequence_hashes=None,
):
    """Return canonical hashes from the chromosome cache or a batch fallback."""
    region_start = max(1, int(region_start))
    region_end = max(region_start, int(region_end))
    if sequence_hashes is not None:
        hash_start = min(region_start, len(sequence_hashes))
        hash_end = min(
            len(sequence_hashes),
            max(hash_start, region_end - k + 1),
        )
        return np.asarray(sequence_hashes[hash_start:hash_end], dtype=np.int32)

    sequence = extract_region(
        fasta_file=fasta_file,
        chr=seq_id,
        region_start=region_start,
        region_end=region_end,
    )
    if not sequence:
        return np.empty(0, dtype=np.int32)
    return canonical_kmer_hashes(sequence, k)


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


def ntr_prism(region, size, kmer=21, verbose=False):
    kmers_low = forward_kmer_hashes(region, kmer)
    kmer_list_low = kmers_low[1:size]
    histo_list_low = calculate_hash_distances(kmer_list_low)
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
                print(
                    "Not a harmonic relationship; using the most populous "
                    f"distance ({most_populous_dist})"
                )
            # A short incidental spacing must not override an unrelated,
            # strongly supported repeat period. Prefer the lowest peak only
            # when it plausibly explains the dominant peak as a harmonic.
            monomer_size = most_populous_dist

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
            if match is not None:
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
    sequence_hashes=None,
    strong_matrix_evidence=False,
):
    candidate_start = max(0, min(int(coordinates[0]), seq_len))
    candidate_end = max(candidate_start, min(int(coordinates[1]), seq_len))
    candidate_boundaries = (candidate_start, candidate_end)

    # Get the sequence of the specified region, subtracted by 1.5 * window size as a potential buffer region.
    interval = math.ceil(window / 2)
    array_seq_size = candidate_end - candidate_start
    if array_seq_size <= 0:
        return None

    # Check if the buffer is too large for the array size
    if array_seq_size < (2 * window):
        window = max(1, math.ceil(array_seq_size / 4))
        interval = math.ceil(window / 2)

    core_start = max(0, candidate_start + window + interval)
    core_end = min(seq_len, candidate_end - window - interval)
    core_seq_size = core_end - core_start

    if core_seq_size <= 0:
        if verbose:
            print(
                f"Removed candidate {candidate_start}-{candidate_end}: "
                "the boundary-refinement core is empty."
            )
        return None

    core_seq = extract_region(
        fasta_file=fasta_file, chr=seq_id, region_start=core_start, region_end=core_end
    )

    if not core_seq:
        if verbose:
            print(
                f"Removed candidate {candidate_start}-{candidate_end}: "
                f"unable to extract core sequence {core_start}-{core_end}."
            )
        return None

    user_prism_res = ntr_prism(core_seq, core_seq_size, k)
    k6_prism_res = (
        user_prism_res if int(k) == 6 else ntr_prism(core_seq, core_seq_size, 6)
    )
    prism_res = (
        user_prism_res
        if user_prism_res[0] not in (None, 0)
        else k6_prism_res
    )
    # Matrix and distal evidence can nominate candidates, but cannot replace
    # direct repeat-period evidence for the candidate itself.
    if prism_res[0] in (None, 0):
        if verbose:
            evidence = "strong" if strong_matrix_evidence else "weak"
            rejected_k = "k=6" if int(k) == 6 else f"k={k} and k=6"
            print(
                f"Removed candidate {candidate_start}-{candidate_end}: "
                f"NTRPrism rejected at {rejected_k}; matrix support was "
                f"{evidence}."
            )
        return None

    if sequence_hashes is None:
        core_hashes = canonical_kmer_hashes(core_seq, k)
    else:
        core_hashes = region_canonical_hashes(
            fasta_file,
            seq_id,
            core_start,
            core_end,
            k,
            sequence_hashes,
        )
    array_kmer_set = set(core_hashes.tolist())
    recurrent_core_kmers = recurrent_kmer_set(core_hashes)

    # This survived the ntr_prism purge, find the boundary
    left_boundary = detect_left_boundary(
        fasta_file=fasta_file,
        array_seq=core_seq,
        array_seq_size=core_seq_size,
        seq_id=seq_id,
        boundary_point=candidate_start,
        k=k,
        window=window,
        interval=interval,
        boundary_chunk_size=100,
        verbosity=verbose,
        prev_border_coordinate=previous_coordinates[1],
        boundary=0,
        expanded=False,
        array_kmer_set=array_kmer_set,
        sequence_hashes=sequence_hashes,
    )

    right_boundary = detect_right_boundary(
        fasta_file=fasta_file,
        array_seq=core_seq,
        array_seq_size=core_seq_size,
        seq_id=seq_id,
        boundary_point=candidate_end,
        limit=None,
        k=k,
        window=window,
        interval=interval,
        boundary_chunk_size=100,
        verbosity=verbose,
        bordering=False,
        boundary=seq_len,
        array_kmer_set=array_kmer_set,
        fallback_array_kmer_set=recurrent_core_kmers,
        sequence_hashes=sequence_hashes,
    )

    if left_boundary is None:
        left_boundary = candidate_start
    if right_boundary is None:
        right_boundary = candidate_end

    left_boundary = max(0, min(int(left_boundary), seq_len))
    right_boundary = max(0, min(int(right_boundary), seq_len))
    if right_boundary <= left_boundary:
        left_boundary, right_boundary = candidate_boundaries

    if classify:
        array_seq_size = right_boundary - left_boundary
        array_hashes = region_canonical_hashes(
            fasta_file,
            seq_id,
            left_boundary,
            right_boundary,
            k,
            sequence_hashes,
        )
        if len(array_hashes) == 0:
            best_match = "Unknown"
        else:
            query_set = set(array_hashes[:array_seq_size].tolist())
            try:
                best_match, results = classify_kmers(query_set, classify, False)
            except Exception:
                best_match = "Unknown"
            if not best_match:
                best_match = "Unknown"

    if verbose:
        evidence = "strong" if strong_matrix_evidence else "weak"
        classification_summary = (
            f" (classification={best_match})" if classify else ""
        )
        print(
            f"Accepted candidate {candidate_start}-{candidate_end}"
            f"{classification_summary}:"
        )
        print(f"    {evidence.capitalize()} matrix support")
        print(
            "    "
            + format_ntrprism_verbose_value(6, k6_prism_res)
        )
        print(
            "    "
            + format_ntrprism_verbose_value(
                k, user_prism_res, user_selected=True
            )
        )
        print(f"    Estimated left boundary: {left_boundary}")
        print(f"    Estimated right boundary: {right_boundary}")

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
    array_kmer_set=None,
    sequence_hashes=None,
):
    boundary_start = boundary_point - window - (interval * 4)
    # if prev_border_coordinate > (boundary_start - interval):
    # print(f"Potential conflict as previous boundary at {prev_border_coordinate} close to {boundary_start}")

    # Ensure the boundary start isn't less than 0
    if boundary_start < boundary:
        boundary_start = boundary

    boundary_end = boundary_point + window + (interval * 4)
    boundary_size = boundary_end - boundary_start
    if boundary_size <= 0:
        return None

    if array_kmer_set is None:
        array_kmer_set = set(canonical_kmer_hashes(array_seq, k).tolist())
    border_kmer_list = region_canonical_hashes(
        fasta_file,
        seq_id,
        boundary_start,
        boundary_end,
        k,
        sequence_hashes,
    ).tolist()
    if not border_kmer_list:
        return None

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
    if last_nonzero_step is None:
        return None
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
                    array_kmer_set=array_kmer_set,
                    sequence_hashes=sequence_hashes,
                )
            else:
                if prev_border_coordinate > boundary_point - window:
                    return prev_border_coordinate + 1
                return max(boundary, boundary_point)

    # Case where boundary needs to be extended to the right
    elif last_nonzero_step >= steps - 1:
        updated_boundary_start = boundary_end - interval
        updated_boundary_end = boundary_end + (interval * 4)
        updated_boundary_start = max(boundary, updated_boundary_start)
        updated_boundary_size = updated_boundary_end - updated_boundary_start
        updated_border_kmer_list = region_canonical_hashes(
            fasta_file,
            seq_id,
            updated_boundary_start,
            updated_boundary_end,
            k,
            sequence_hashes,
        ).tolist()
        if not updated_border_kmer_list:
            return None
        updated_steps = math.ceil(updated_boundary_size / 100)
        last_nonzero_step = find_target_fixed_window_left(
            array_kmer_set=array_kmer_set,
            border_kmer_list=updated_border_kmer_list,
            boundary_size=updated_boundary_size,
            offset=updated_boundary_start,
            verbose=verbosity,
            step_size=updated_steps,
        )
        if last_nonzero_step is None:
            return None

        # Error case
        if last_nonzero_step <= 0 or last_nonzero_step >= updated_steps - 1:
            return None

        # find_target_fixed_window_left(array_kmer_set, updated_border_kmer_list, updated_boundary_size, updated_boundary_start, updated_steps)
        possible_range = (
            updated_boundary_start + (100 * last_nonzero_step),
            updated_boundary_start + (100 * last_nonzero_step) + 200,
        )
        # print(f"New possible range: {possible_range[0]}-{possible_range[1]}")
        border_array_start = boundary_chunk_size * last_nonzero_step
        border_array_end = border_array_start + (2 * boundary_chunk_size)

        estimated_index = find_last_matching_index(
            updated_border_kmer_list,
            array_kmer_set,
            border_array_start,
            border_array_end,
        )
        return max(boundary, estimated_index + updated_boundary_start + k)

    else:
        border_array_start = boundary_chunk_size * last_nonzero_step
        border_array_end = border_array_start + (2 * boundary_chunk_size)

        estimated_index = find_last_matching_index(
            border_kmer_list, array_kmer_set, border_array_start, border_array_end
        )
        return max(
            boundary, estimated_index + boundary_start - k
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
    array_kmer_set=None,
    fallback_array_kmer_set=None,
    sequence_hashes=None,
):
    boundary_start = max(0, boundary_point - window - interval)
    boundary_end = boundary_point + window + interval
    # Keep the search inside the sequence.
    if boundary_end > boundary:
        boundary_end = boundary
    boundary_size = boundary_end - boundary_start
    if boundary_size <= 0:
        return None
    if array_kmer_set is None:
        array_kmer_set = set(canonical_kmer_hashes(array_seq, k).tolist())
    border_kmer_list = region_canonical_hashes(
        fasta_file,
        seq_id,
        boundary_start,
        boundary_end,
        k,
        sequence_hashes,
    ).tolist()
    if not border_kmer_list:
        return None

    def search_left_of_initial_range(search_kmers):
        """Recover a boundary that lies before the initial search window."""
        if not search_kmers:
            return None
        updated_boundary_end = min(boundary, boundary_start + interval)
        updated_boundary_start = max(0, boundary_start - (interval * 4))
        updated_boundary_size = updated_boundary_end - updated_boundary_start
        if updated_boundary_size <= 0:
            return None
        updated_border_kmer_list = region_canonical_hashes(
            fasta_file,
            seq_id,
            updated_boundary_start,
            updated_boundary_end,
            k,
            sequence_hashes,
        ).tolist()
        if not updated_border_kmer_list:
            return None
        updated_steps = math.ceil(updated_boundary_size / boundary_chunk_size)
        updated_last_nonzero_step = find_target_fixed_window_right(
            array_kmer_set=search_kmers,
            border_kmer_list=updated_border_kmer_list,
            boundary_size=updated_boundary_size,
            offset=updated_boundary_start,
            verbose=verbosity,
            step_size=updated_steps,
        )
        if updated_last_nonzero_step is None:
            return None
        border_array_start = boundary_chunk_size * updated_last_nonzero_step
        border_array_end = border_array_start + (2 * boundary_chunk_size)
        estimated_index = find_last_matching_index(
            updated_border_kmer_list,
            search_kmers,
            border_array_start,
            border_array_end,
        )
        estimate = estimated_index + updated_boundary_start + k
        return min(boundary, max(0, estimate))

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
    if last_nonzero_step is None:
        return search_left_of_initial_range(fallback_array_kmer_set)
    possible_range = (
        boundary_start + (boundary_chunk_size * last_nonzero_step),
        boundary_start
        + (boundary_chunk_size * last_nonzero_step)
        + 2 * boundary_chunk_size,
    )

    # Case for boundary extending to the right
    if last_nonzero_step >= steps - 1:
        if boundary_end >= boundary:
            return boundary

        updated_boundary_start = max(0, boundary_end - interval)
        updated_boundary_end = min(boundary, boundary_end + (interval * 4))
        updated_boundary_size = updated_boundary_end - updated_boundary_start
        if updated_boundary_size <= 0:
            return None
        updated_border_kmer_list = region_canonical_hashes(
            fasta_file,
            seq_id,
            updated_boundary_start,
            updated_boundary_end,
            k,
            sequence_hashes,
        ).tolist()
        if not updated_border_kmer_list:
            return None
        updated_steps = math.ceil(updated_boundary_size / boundary_chunk_size)
        updated_last_nonzero_step = find_target_fixed_window_right(
            array_kmer_set=array_kmer_set,
            border_kmer_list=updated_border_kmer_list,
            boundary_size=updated_boundary_size,
            offset=updated_boundary_start,
            verbose=verbosity,
            step_size=updated_steps,
        )
        if updated_last_nonzero_step is None:
            return None
        if (
            updated_last_nonzero_step >= updated_steps - 1
            and updated_boundary_end < boundary
        ):
            return None

        border_array_start = boundary_chunk_size * updated_last_nonzero_step
        border_array_end = border_array_start + (2 * boundary_chunk_size)
        estimated_index = find_last_matching_index(
            updated_border_kmer_list,
            array_kmer_set,
            border_array_start,
            border_array_end,
        )
        return min(boundary, estimated_index + updated_boundary_start + k)

    # Case for boundary extending left
    elif last_nonzero_step <= 1:
        return search_left_of_initial_range(array_kmer_set)

    else:
        border_array_start = boundary_chunk_size * last_nonzero_step
        border_array_end = border_array_start + (2 * boundary_chunk_size)

        estimated_index = find_last_matching_index(
            border_kmer_list, array_kmer_set, border_array_start, border_array_end
        )
        return min(
            boundary, max(0, estimated_index + boundary_start + k)
        )  # Return the estimated boundary position adjusted by k-mer size


def find_target_fixed_window_left(
    array_kmer_set, border_kmer_list, boundary_size, offset, verbose, step_size=100
):
    chunk_count = max(1, int(step_size))
    chunk_size = max(1, math.ceil(boundary_size / chunk_count))
    threshold = math.ceil(chunk_size / 4)
    last_nonzero_step = None

    shared_counts = []

    count = 0
    for i in range(chunk_count - 1, -1, -1):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, len(border_kmer_list))
        border_chunk = border_kmer_list[start_idx:end_idx]
        shared_count = sum(kmer in array_kmer_set for kmer in border_chunk)
        shared_counts.append((i, shared_count))

        if shared_count > threshold and count < 10:
            last_nonzero_step = i
            count = 0
        else:
            count += 1

    plot_height = 10  # Reduced from 20 to 10 for half the height

    show_histogram = (
        verbose
        and last_nonzero_step is not None
        and 1 < last_nonzero_step < chunk_count - 1
    )
    if show_histogram:
        graph_columns, max_count = boundary_graph_columns(
            shared_counts, last_nonzero_step, plot_height
        )
        for y in range(plot_height, -1, -1):
            line = f"{str(int(y * max_count / plot_height)).rjust(4)} | "
            for contains_stop, height in reversed(graph_columns):
                if height >= y:
                    line += "|" if contains_stop else "#"
                else:
                    line += " "
            print(line)

        # X-axis line
        graph_width = len(graph_columns)
        axis_line = "     " + "-" * graph_width
        print(axis_line)

        # X-axis labels: leftmost and rightmost only
        print(boundary_graph_label_line(offset, boundary_size, graph_width))

    return last_nonzero_step


def find_target_fixed_window_right(
    array_kmer_set, border_kmer_list, boundary_size, offset, verbose, step_size=100
):
    chunk_count = max(1, int(step_size))
    chunk_size = max(1, math.ceil(boundary_size / chunk_count))
    threshold = math.ceil(chunk_size / 4)
    last_nonzero_step = None

    shared_counts = []
    counter = 0
    for i in range(chunk_count):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, len(border_kmer_list))
        border_chunk = border_kmer_list[start_idx:end_idx]
        shared_count = sum(kmer in array_kmer_set for kmer in border_chunk)
        shared_counts.append((i, shared_count))

        if shared_count > threshold and counter < 10:
            last_nonzero_step = i
            counter = 0
        else:
            counter += 1

    plot_height = 10  # match left-hand version

    show_histogram = (
        verbose
        and last_nonzero_step is not None
        and 1 < last_nonzero_step < chunk_count - 1
    )
    if show_histogram:
        graph_columns, max_count = boundary_graph_columns(
            shared_counts, last_nonzero_step, plot_height
        )
        for y in range(plot_height, -1, -1):
            line = f"{str(int(y * max_count / plot_height)).rjust(4)} | "
            for contains_stop, height in graph_columns:
                if height >= y:
                    line += "|" if contains_stop else "#"
                else:
                    line += " "
            print(line)

        # X-axis line
        graph_width = len(graph_columns)
        axis_line = "     " + "-" * graph_width
        print(axis_line)

        # X-axis labels: leftmost and rightmost
        print(boundary_graph_label_line(offset, boundary_size, graph_width))

    return last_nonzero_step


def find_last_matching_index(
    border_kmer_list, array_kmer_set, border_array_start, border_array_end
):
    start = max(0, min(int(border_array_start), len(border_kmer_list)))
    end = max(start, min(int(border_array_end), len(border_kmer_list)))
    for i in range(end - 1, start - 1, -1):
        if border_kmer_list[i] in array_kmer_set:
            return i
    return start


def candidates_share_band_border(current_end, next_start, band, tolerance):
    """Return True only when adjacent candidates touch the same band edge."""
    if band <= 0:
        return False
    boundary = round(current_end / band) * band
    if boundary <= 0:
        return False
    return (
        abs(current_end - boundary) <= tolerance
        and abs(next_start - boundary) <= tolerance
        and next_start >= current_end - tolerance
    )


def report_borders(
    fa,
    seq_id: str,
    seq_len: int,
    band: int,
    offset: int,
    window: int,
    k: int,
    df: pl.DataFrame,
    classify: "str | bool",
    verbose: bool,
    quiet: bool,
    sequence_hashes=None,
    candidate_windows=None,
    candidate_matrix_support=None,
) -> None:
    """
    Infer satellite locations & precise boundaries for regions described in `df`.
    """
    # Convert band (in millions) to base-pairs without truncating sub-Mbp bands.
    band = int(float(band) * 1_000_000)
    if band <= 0:
        raise ValueError("band must be greater than zero")

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
        displayed_windows = sorted(
            set([window] if candidate_windows is None else candidate_windows)
        )
        print(
            f"Processing {seq_id} boundaries with window size(s) "
            f"{', '.join(map(str, displayed_windows))} and k-mer size {k}:\n"
        )

    # Prepare starts/ends relative to the provided offset
    starts = [s - offset for s in df["start"].to_list()]
    ends = [e - offset for e in df["end"].to_list()]

    # Ensure non-negative
    starts = [s if s > 0 else 0 for s in starts]
    ends = [e if e > 0 else 0 for e in ends]
    if candidate_windows is None:
        refinement_windows = [window] * len(starts)
    else:
        refinement_windows = [int(value) for value in candidate_windows]
        if len(refinement_windows) != len(starts):
            raise ValueError(
                "candidate_windows must contain one value for every candidate"
            )
        if any(value <= 0 for value in refinement_windows):
            raise ValueError("candidate window sizes must be positive")
    if candidate_matrix_support is None:
        strong_matrix_evidence = [False] * len(starts)
    else:
        strong_matrix_evidence = [bool(value) for value in candidate_matrix_support]
        if len(strong_matrix_evidence) != len(starts):
            raise ValueError(
                "candidate_matrix_support must contain one value for every candidate"
            )

    new_starts = []
    new_ends = []
    classification = []
    monomer = []
    periodicity = []
    hor = []
    previous_new_boundaries = []

    assert len(starts) == len(ends)

    def append_result(result):
        start, end, name, prism_result = result
        if prism_result is None or len(prism_result) < 4:
            raise ValueError("Boundary result is missing repeat metadata")
        start = int(start)
        end = int(end)
        monomer_value = prism_result[0]
        if monomer_value in (None, 0):
            raise ValueError("Boundary result has no valid NTRPrism monomer")
        periodicity_value = prism_result[2]
        hor_value = bool(prism_result[1])
        new_starts.append(start)
        new_ends.append(end)
        classification.append(name)
        monomer.append(monomer_value)
        periodicity.append(periodicity_value)
        hor.append(hor_value)
        previous_new_boundaries.append((start, end))

    i = 0
    while i < len(starts):
        next_i = i + 1
        previous_coordinates = (
            previous_new_boundaries[-1] if previous_new_boundaries else (0, 1)
        )
        coordinates = (int(starts[i]) + 1, int(ends[i]) + 1)
        refinement_window = refinement_windows[i]
        candidate_has_strong_matrix_evidence = strong_matrix_evidence[i]

        if next_i < len(starts) and candidates_share_band_border(
            current_end=ends[i],
            next_start=starts[next_i],
            band=band,
            tolerance=max(refinement_window, refinement_windows[next_i]) * 3,
        ):
            coordinates = (coordinates[0], int(ends[next_i]) + 1)
            refinement_window = min(
                refinement_window, refinement_windows[next_i]
            )
            candidate_has_strong_matrix_evidence = (
                candidate_has_strong_matrix_evidence
                or strong_matrix_evidence[next_i]
            )
            del starts[next_i]
            del ends[next_i]
            del refinement_windows[next_i]
            del strong_matrix_evidence[next_i]

        if verbose:
            print("\n" + "=" * BOUNDARY_GRAPH_MAX_WIDTH)
            print(
                "Analyzing candidate: "
                f"left boundary={coordinates[0]}, right boundary={coordinates[1]}"
            )

        try:
            updated_boundaries = detect_precise_boundaries(
                fasta_file=fa,
                seq_id=seq_id,
                seq_len=seq_len,
                window=refinement_window,
                k=k,
                coordinates=coordinates,
                verbose=verbose,
                classify=supersets_dict,
                previous_coordinates=previous_coordinates,
                sequence_hashes=sequence_hashes,
                strong_matrix_evidence=candidate_has_strong_matrix_evidence,
            )
            # A None result is an intentional ntr_prism rejection, not an
            # extension error, so it remains excluded.
            if updated_boundaries is None:
                i += 1
                continue
            append_result(updated_boundaries)
        except Exception as error:
            if not quiet:
                print(
                    f"Removed candidate {coordinates[0]}-{coordinates[1]}: "
                    f"boundary refinement failed ({error})."
                )

        i += 1

    return new_starts, new_ends, classification, monomer, periodicity, hor
