from numba import njit, prange
import numpy as np


@njit
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


@njit(parallel=True)
def intersection_matrix(overlapping, non_overlapping, k):
    n = len(overlapping)
    mat = np.empty((n, n), dtype=np.float64)
    powk = 1.0 / k
    for i in prange(n):
        a = non_overlapping[i]
        a_prime = overlapping[i]
        len_a = len(a)
        inv_len_a = 1.0 / len_a
        for j in range(i, n):
            b = non_overlapping[j]
            b_prime = overlapping[j]
            len_b = len(b)
            if len_a == 0 or len_b == 0:
                mat[i, j] = 0.0
                mat[j, i] = 0.0
                continue
            inv_len_b = 1.0 / len_b
            inter1 = intersection_len(a, b_prime) * inv_len_a
            inter2 = intersection_len(a_prime, b) * inv_len_b
            value = max(inter1, inter2)
            score = (value**powk) * 100.0
            mat[i, j] = score
            mat[j, i] = score
    return mat


@njit(parallel=True)
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
