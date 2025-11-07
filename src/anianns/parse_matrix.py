from collections import Counter
import numpy as np
from scipy import ndimage

def append_coordinates(satellite_coordinate_list, prefix, matrix, window, off_diagonal):
    M_diag, M_offdiag = split_diagonal_attached(matrix)
    if off_diagonal:
        merged = merge_shared_boundaries(M_offdiag, prefix, False)
    else:
        merged = merge_shared_boundaries(matrix, prefix, False)
        for coordinates in merged:
            satellite_coordinate_list.append(coordinates)
    return satellite_coordinate_list


def check_same_start_end(pairs, s, start, window=2000):
    lo, hi = s - window, s + window
    if start:
        return [(x,y) for x,y in pairs if lo <= x <= hi]
    else:
        return [(x,y) for x,y in pairs if lo <= y <= hi]

def check_contained(pairs, start, end):
    return [(x,y) for x,y in pairs if (x < start) and (y > end)]


def get_diagonal_span(matrix, window, zero_tol):
    n = matrix.shape[0]
    lengths = [0] * n
    coords  = [(0, 0)] * n

    for i in range(n):
        row = matrix[i]         # local view of row i
        start = end = i
        cnt = 0

        # scan left of diagonal
        zeros = 0
        for j in range(i - 1, -1, -1):
            if row[j] != 0:
                cnt += 1
                start = j
                zeros = 0
            else:
                zeros += 1
                if zeros == zero_tol:
                    break

        # scan right of diagonal
        zeros = 0
        for j in range(i + 1, n):
            if row[j] != 0:
                cnt += 1
                end = j
                zeros = 0
            else:
                zeros += 1
                if zeros == zero_tol:
                    break

        lengths[i] = cnt
        coords[i]  = (start * window, (end * window) + window)
    
    tuple_counts = Counter(coords)
    sorted_items = sorted(
        ((k, v) for k, v in tuple_counts.items() if v >= 3),
        key=lambda item: item[0][0]
    )
    return sorted_items

def merge_shared_boundaries(intervals, prefix, verbose=False):
    """
    intervals: list of ((start, end), count)
    returns: list of (start, end) with “one‑off” or containment conflicts removed
    """
    out = []
    # cache the conflict‑check functions for faster lookup
    chk_same = check_same_start_end
    chk_cont  = check_contained

    for (x, y), _count in intervals:
        c_x = chk_same(out, x, True)
        c_y = chk_same(out, y, False)
        c_c = chk_cont(out, x, y)
        has_conflict = c_x or c_y or c_c

        if verbose:
            if has_conflict:
                reasons = []
                if c_x: reasons.append("one‑off start")
                if c_y: reasons.append("one‑off end")
                if c_c: reasons.append("contained interval")
                #print(f"Conflict(s) for ({x}, {y}): {', '.join(reasons)}")
            else:
                #print(f"No conflicts for ({x}, {y}) – keeping it.")
                out.append((x, y))
        else:
            if not has_conflict:
                out.append((x, y))
    return [(x + prefix, y + prefix) for x, y in out]

def split_diagonal_attached(M):
    binary = M > 0
    vertical_structure = np.array([[0, 1, 0],
                          [0, 1, 0],
                          [0, 1, 0]], dtype=int)
    horizontal_structure = np.array([[0, 0, 0],
                          [1, 1, 1],
                          [0, 0, 0]], dtype=int)
    vertical_labeled, _ = ndimage.label(binary, structure=vertical_structure)
    horizontal_labeled, _ = ndimage.label(binary, structure=horizontal_structure)

    diag_indices = np.arange(min(M.shape))
    vertical_diag_labels = np.unique(vertical_labeled[diag_indices, diag_indices])
    horizontal_diag_labels = np.unique(horizontal_labeled[diag_indices, diag_indices])
    
    vertical_keep_mask = np.isin(vertical_labeled, vertical_diag_labels)
    horizontal_keep_mask = np.isin(horizontal_labeled, horizontal_diag_labels)

    combined_mask = np.minimum(vertical_keep_mask, horizontal_keep_mask)

    M_diag = M * combined_mask
    M_offdiag = M * (~combined_mask)

    return M_diag, M_offdiag
