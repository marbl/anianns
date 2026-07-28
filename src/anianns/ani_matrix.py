from numba import njit, prange
import numpy as np


@njit(cache=True)
def intersection_len(a, b):
    """Two-pointer intersection count of sorted int arrays."""
    count = 0
    i = j = 0
    len_a = len(a)
    len_b = len(b)
    while i < len_a and j < len_b:
        if a[i] == b[j]:
            count += 1
            i += 1
            j += 1
        elif a[i] < b[j]:
            i += 1
        else:
            j += 1
    return count


@njit(cache=True, inline="always")
def intersection_reaches_threshold(a, b, required_count):
    """Return early once an intersection must pass or can no longer pass."""
    count = 0
    i = j = 0
    len_a = len(a)
    len_b = len(b)
    while i < len_a and j < len_b:
        if a[i] == b[j]:
            count += 1
            i += 1
            j += 1
            if count >= required_count:
                return True
        elif a[i] < b[j]:
            i += 1
        else:
            j += 1

        # Even matching every remaining value cannot reach the threshold.
        if count + min(len_a - i, len_b - j) < required_count:
            return False
    return count >= required_count


@njit(cache=True, inline="always")
def passes_identity(a, a_prime, b, b_prime, minimum_similarity):
    """Test one symmetric ANI cell without calculating or storing its score."""
    len_a = len(a)
    len_b = len(b)
    if len_a == 0 or len_b == 0:
        return False
    if intersection_reaches_threshold(a, b_prime, minimum_similarity * len_a):
        return True
    return intersection_reaches_threshold(
        a_prime, b, minimum_similarity * len_b
    )


@njit(cache=True)
def _triangular_block_coordinates(n, block_size):
    """Build similarly sized upper-triangle work units for parallel kernels."""
    block_count = (n + block_size - 1) // block_size
    pair_count = block_count * (block_count + 1) // 2
    block_rows = np.empty(pair_count, dtype=np.int64)
    block_columns = np.empty(pair_count, dtype=np.int64)
    pair_index = 0
    for block_row in range(block_count):
        for block_column in range(block_row, block_count):
            block_rows[pair_index] = block_row
            block_columns[pair_index] = block_column
            pair_index += 1
    return block_rows, block_columns


@njit(cache=True, parallel=True)
def intersection_matrix(overlapping, non_overlapping, k):
    n = len(overlapping)
    mat = np.empty((n, n), dtype=np.float64)
    powk = 1.0 / k
    block_size = 32
    block_rows, block_columns = _triangular_block_coordinates(n, block_size)
    for block_index in prange(len(block_rows)):
        row_start = block_rows[block_index] * block_size
        row_stop = min(row_start + block_size, n)
        column_start = block_columns[block_index] * block_size
        column_stop = min(column_start + block_size, n)
        for i in range(row_start, row_stop):
            a = non_overlapping[i]
            a_prime = overlapping[i]
            len_a = len(a)
            first_column = max(i, column_start)
            for j in range(first_column, column_stop):
                b = non_overlapping[j]
                b_prime = overlapping[j]
                len_b = len(b)
                if len_a == 0 or len_b == 0:
                    mat[i, j] = 0.0
                    mat[j, i] = 0.0
                    continue
                inv_len_a = 1.0 / len_a
                inv_len_b = 1.0 / len_b
                inter1 = intersection_len(a, b_prime) * inv_len_a
                inter2 = intersection_len(a_prime, b) * inv_len_b
                value = max(inter1, inter2)
                score = (value**powk) * 100.0
                mat[i, j] = score
                mat[j, i] = score
    return mat


@njit(cache=True, parallel=True)
def intersection_matrix_with_threshold(overlapping, non_overlapping, k, identity):
    """Return exact identities and the threshold predicate in one traversal."""
    n = len(overlapping)
    identities = np.empty((n, n), dtype=np.float64)
    thresholded = np.zeros((n, n), dtype=np.bool_)
    powk = 1.0 / k
    minimum_similarity = (identity / 100.0) ** k
    block_size = 32
    block_rows, block_columns = _triangular_block_coordinates(n, block_size)
    for block_index in prange(len(block_rows)):
        row_start = block_rows[block_index] * block_size
        row_stop = min(row_start + block_size, n)
        column_start = block_columns[block_index] * block_size
        column_stop = min(column_start + block_size, n)
        for i in range(row_start, row_stop):
            a = non_overlapping[i]
            a_prime = overlapping[i]
            len_a = len(a)
            first_column = max(i, column_start)
            for j in range(first_column, column_stop):
                b = non_overlapping[j]
                b_prime = overlapping[j]
                len_b = len(b)
                if len_a == 0 or len_b == 0:
                    identities[i, j] = 0.0
                    identities[j, i] = 0.0
                    continue
                inv_len_a = 1.0 / len_a
                inv_len_b = 1.0 / len_b
                intersection1 = intersection_len(a, b_prime)
                intersection2 = intersection_len(a_prime, b)
                inter1 = intersection1 * inv_len_a
                inter2 = intersection2 * inv_len_b
                similarity = max(inter1, inter2)
                score = (similarity**powk) * 100.0
                identities[i, j] = score
                identities[j, i] = score
                if (
                    intersection1 >= minimum_similarity * len_a
                    or intersection2 >= minimum_similarity * len_b
                ):
                    thresholded[i, j] = True
                    thresholded[j, i] = True
    return identities, thresholded


@njit(cache=True, parallel=True)
def intersection_matrix_rectangular(
    left_overlapping,
    left_non_overlapping,
    right_overlapping,
    right_non_overlapping,
    k,
):
    """Calculate exact ANI scores between two different window groups."""
    n_left = len(left_overlapping)
    n_right = len(right_overlapping)
    mat = np.empty((n_left, n_right), dtype=np.float64)
    powk = 1.0 / k
    for i in prange(n_left):
        a = left_non_overlapping[i]
        a_prime = left_overlapping[i]
        len_a = len(a)
        for j in range(n_right):
            b = right_non_overlapping[j]
            b_prime = right_overlapping[j]
            len_b = len(b)
            if len_a == 0 or len_b == 0:
                mat[i, j] = 0.0
                continue
            inter1 = intersection_len(a, b_prime) / len_a
            inter2 = intersection_len(a_prime, b) / len_b
            mat[i, j] = (max(inter1, inter2) ** powk) * 100.0
    return mat


@njit(cache=True, parallel=True)
def intersection_matrix_thresholded(overlapping, non_overlapping, k, identity):
    """Return a compact 0/1 matrix containing only scores above ``identity``."""
    n = len(overlapping)
    mat = np.zeros((n, n), dtype=np.bool_)
    minimum_similarity = (identity / 100.0) ** k
    block_size = 32
    block_rows, block_columns = _triangular_block_coordinates(n, block_size)
    for block_index in prange(len(block_rows)):
        row_start = block_rows[block_index] * block_size
        row_stop = min(row_start + block_size, n)
        column_start = block_columns[block_index] * block_size
        column_stop = min(column_start + block_size, n)
        for i in range(row_start, row_stop):
            a = non_overlapping[i]
            a_prime = overlapping[i]
            first_column = max(i, column_start)
            for j in range(first_column, column_stop):
                b = non_overlapping[j]
                b_prime = overlapping[j]
                if passes_identity(a, a_prime, b, b_prime, minimum_similarity):
                    mat[i, j] = True
                    mat[j, i] = True
    return mat


@njit(cache=True, parallel=True)
def intersection_matrix_cross_groups_thresholded(
    overlapping,
    non_overlapping,
    group_ids,
    k,
    identity,
    skip_adjacent_groups=False,
):
    """Compare only windows from different precomputed band groups."""
    n = len(overlapping)
    mat = np.zeros((n, n), dtype=np.bool_)
    minimum_similarity = (identity / 100.0) ** k
    block_size = 32
    block_rows, block_columns = _triangular_block_coordinates(n, block_size)
    for block_index in prange(len(block_rows)):
        row_start = block_rows[block_index] * block_size
        row_stop = min(row_start + block_size, n)
        column_start = block_columns[block_index] * block_size
        column_stop = min(column_start + block_size, n)
        for i in range(row_start, row_stop):
            a = non_overlapping[i]
            a_prime = overlapping[i]
            first_column = max(i + 1, column_start)
            for j in range(first_column, column_stop):
                if group_ids[i] == group_ids[j]:
                    continue
                if (
                    skip_adjacent_groups
                    and abs(group_ids[i] - group_ids[j]) == 1
                ):
                    continue
                if passes_identity(
                    a,
                    a_prime,
                    non_overlapping[j],
                    overlapping[j],
                    minimum_similarity,
                ):
                    mat[i, j] = True
                    mat[j, i] = True
    return mat


@njit(cache=True, parallel=True)
def intersection_matrix_rectangular_thresholded(
    row_overlapping,
    row_non_overlapping,
    row_indices,
    column_overlapping,
    column_non_overlapping,
    column_indices,
    k,
    identity,
):
    """Compare selected windows from two different sequence bands."""
    row_count = len(row_indices)
    column_count = len(column_indices)
    matrix = np.zeros((row_count, column_count), dtype=np.bool_)
    minimum_similarity = (identity / 100.0) ** k
    for selected_row in prange(row_count):
        row_index = row_indices[selected_row]
        a = row_non_overlapping[row_index]
        a_prime = row_overlapping[row_index]
        for selected_column in range(column_count):
            column_index = column_indices[selected_column]
            if passes_identity(
                a,
                a_prime,
                column_non_overlapping[column_index],
                column_overlapping[column_index],
                minimum_similarity,
            ):
                matrix[selected_row, selected_column] = True
    return matrix


@njit(cache=True)
def _diagonal_row_bounds(
    row_index,
    overlapping,
    non_overlapping,
    minimum_similarity,
    zero_tolerance,
):
    a = non_overlapping[row_index]
    a_prime = overlapping[row_index]
    left_bound = row_index
    right_bound = row_index

    zeros = 0
    for column in range(row_index - 1, -1, -1):
        if passes_identity(
            a,
            a_prime,
            non_overlapping[column],
            overlapping[column],
            minimum_similarity,
        ):
            left_bound = column
            zeros = 0
        else:
            zeros += 1
            if zeros == zero_tolerance:
                break

    zeros = 0
    for column in range(row_index + 1, len(overlapping)):
        if passes_identity(
            a,
            a_prime,
            non_overlapping[column],
            overlapping[column],
            minimum_similarity,
        ):
            right_bound = column
            zeros = 0
        else:
            zeros += 1
            if zeros == zero_tolerance:
                break

    return left_bound, right_bound


@njit(cache=True, parallel=True)
def diagonal_span_bounds(
    overlapping, non_overlapping, k, identity, zero_tolerance
):
    """Find each diagonal row span without materializing the full matrix."""
    n = len(overlapping)
    starts = np.empty(n, dtype=np.int64)
    ends = np.empty(n, dtype=np.int64)
    minimum_similarity = (identity / 100.0) ** k

    for row_index in prange(n):
        starts[row_index], ends[row_index] = _diagonal_row_bounds(
            row_index,
            overlapping,
            non_overlapping,
            minimum_similarity,
            zero_tolerance,
        )

    return starts, ends


@njit(cache=True, parallel=True)
def periodic_lag_matches(
    overlapping,
    non_overlapping,
    k,
    identity,
    max_lag,
):
    """Compare only bounded diagonal offsets for periodic-stripe detection."""
    n = len(overlapping)
    bounded_lag = min(max(0, max_lag), max(0, n - 1))
    matches = np.zeros((bounded_lag + 1, n), dtype=np.bool_)
    minimum_similarity = (identity / 100.0) ** k
    for lag in prange(1, bounded_lag + 1):
        for row_index in range(n - lag):
            column_index = row_index + lag
            if passes_identity(
                non_overlapping[row_index],
                overlapping[row_index],
                non_overlapping[column_index],
                overlapping[column_index],
                minimum_similarity,
            ):
                matches[lag, row_index] = True
    return matches


@njit(cache=True, parallel=True)
def intersection_matrix_selected_thresholded(
    overlapping, non_overlapping, selected_indices, k, identity
):
    """Materialize only selected candidate-neighborhood windows."""
    selected_count = len(selected_indices)
    matrix = np.zeros((selected_count, selected_count), dtype=np.bool_)
    minimum_similarity = (identity / 100.0) ** k
    for selected_i in prange(selected_count):
        i = selected_indices[selected_i]
        a = non_overlapping[i]
        a_prime = overlapping[i]
        for selected_j in range(selected_i, selected_count):
            j = selected_indices[selected_j]
            if passes_identity(
                a,
                a_prime,
                non_overlapping[j],
                overlapping[j],
                minimum_similarity,
            ):
                matrix[selected_i, selected_j] = True
                matrix[selected_j, selected_i] = True
    return matrix


@njit(cache=True, parallel=True)
def intersection_matrix_selected_vs_all_thresholded(
    overlapping, non_overlapping, selected_indices, k, identity
):
    """Compare selected candidate rows against every window in one band."""
    selected_count = len(selected_indices)
    window_count = len(overlapping)
    matrix = np.zeros((selected_count, window_count), dtype=np.bool_)
    minimum_similarity = (identity / 100.0) ** k
    for selected_row in prange(selected_count):
        row_index = selected_indices[selected_row]
        a = non_overlapping[row_index]
        a_prime = overlapping[row_index]
        for column in range(window_count):
            if passes_identity(
                a,
                a_prime,
                non_overlapping[column],
                overlapping[column],
                minimum_similarity,
            ):
                matrix[selected_row, column] = True
    return matrix


@njit(cache=True, parallel=True)
def intersection_matrix_inverted(
    A, B, overlapping_A, non_overlapping_A, overlapping_B, non_overlapping_B, k
):
    # Assume both matrices are square
    merged_size = A.shape[0] + B.shape[0]
    merged_matrix = np.zeros((merged_size, merged_size), dtype=A.dtype)
    # Place A in top-left
    merged_matrix[: A.shape[0], : A.shape[1]] = A
    # Place B in bottom-right
    merged_matrix[A.shape[0] :, A.shape[1] :] = B
    # Iterate over everything
    for i in prange(A.shape[0], merged_size):
        set_b = non_overlapping_B[i - A.shape[0]]
        set_b_prime = overlapping_B[i - A.shape[0]]
        len_b = len(set_b)
        for j in range((i - B.shape[0]), A.shape[0]):
            set_a = non_overlapping_A[j - A.shape[0]]
            set_a_prime = overlapping_A[j - A.shape[0]]
            len_a = len(set_a)
            if len_a == 0 or len_b == 0:
                merged_matrix[i, j] = 0.0
                merged_matrix[j, i] = 0.0
                continue
            inter1 = intersection_len(set_a, set_b_prime) / len_a
            inter2 = intersection_len(set_a_prime, set_b) / len_b
            value = max(inter1, inter2)
            merged_matrix[i, j] = value ** (1 / k) * 100
            merged_matrix[j, i] = value ** (1 / k) * 100
    return merged_matrix
