import numpy as np

from anianns.ani_matrix import (
    intersection_len,
    intersection_matrix,
    intersection_matrix_inverted,
)


def test_intersection_len_counts_sorted_overlaps():
    assert intersection_len.py_func(np.array([1, 2, 4]), np.array([2, 3, 4])) == 2


def test_intersection_matrix_is_symmetric_and_handles_empty_windows():
    overlapping = [
        np.array([1, 2, 3], dtype=np.int32),
        np.array([2, 3, 4], dtype=np.int32),
        np.array([4, 5], dtype=np.int32),
    ]
    non_overlapping = [
        np.array([1, 2], dtype=np.int32),
        np.array([2, 4], dtype=np.int32),
        np.array([4, 5], dtype=np.int32),
    ]

    matrix = intersection_matrix.py_func(overlapping, non_overlapping, 2)

    assert matrix.shape == (3, 3)
    assert np.allclose(matrix, matrix.T)
    assert matrix[0, 0] == 100.0
    assert matrix[2, 2] == 100.0
    assert matrix[0, 1] > 0.0


def test_intersection_matrix_inverted_merges_two_blocks():
    a = np.array([[100.0, 80.0], [80.0, 100.0]])
    b = np.array([[100.0, 75.0], [75.0, 100.0]])
    overlapping_a = [np.array([1, 2], dtype=np.int32), np.array([2, 3], dtype=np.int32)]
    non_overlapping_a = [np.array([1, 2], dtype=np.int32), np.array([2, 3], dtype=np.int32)]
    overlapping_b = [np.array([2, 4], dtype=np.int32), np.array([3, 4], dtype=np.int32)]
    non_overlapping_b = [np.array([2, 4], dtype=np.int32), np.array([3, 4], dtype=np.int32)]

    merged = intersection_matrix_inverted.py_func(
        a,
        b,
        overlapping_a,
        non_overlapping_a,
        overlapping_b,
        non_overlapping_b,
        2,
    )

    assert merged.shape == (4, 4)
    assert np.allclose(np.diag(merged), 100.0)
    assert np.allclose(merged, merged.T)
    assert merged[2, 0] > 0.0
