from collections import Counter
import numpy as np
from scipy import ndimage
from statistics import median
import matplotlib.pyplot as plt


def append_coordinates(satellite_coordinate_list, prefix, matrix, window, off_diagonal):
    M_diag, M_offdiag = split_diagonal_attached(matrix)
    if off_diagonal:
        merged = merge_shared_boundaries(M_offdiag, prefix, False)
    else:
        merged = merge_shared_boundaries(matrix, prefix, False)
        for coordinates in merged:
            satellite_coordinate_list.append(coordinates)
    return satellite_coordinate_list


def check_same_start_end(pairs, s, window, start):
    lo, hi = s - window, s + window
    if start:
        return [(x, y, count) for x, y, count in pairs if lo <= x <= hi]
    else:
        return [(x, y, count) for x, y, count in pairs if lo <= y <= hi]


def check_new_contained(pairs, start, end):
    return [(x, y, count) for x, y, count in pairs if (x > start) and (y < end)]


def check_new_spans(pairs, start, end):
    return [(x, y, count) for x, y, count in pairs if (x < start) and (y > end)]


def find_in_range_x(data, min_val, max_val):
    return [(i, t) for i, t in enumerate(data) if min_val <= t[0] <= max_val]


def find_in_range_y(data, min_val, max_val):
    return [(i, t) for i, t in enumerate(data) if min_val <= t[1] <= max_val]


def find_contained(data, min_val_x, max_val_x, min_val_y, max_val_y):
    return [
        (i, t) for i, t in enumerate(data) if (min_val_x > t[0]) and (max_val_y < t[1])
    ]


def find_spanning(data, min_val_x, max_val_x, min_val_y, max_val_y):
    return [
        (i, t)
        for i, t in enumerate(data)
        if ((min_val_x < t[0]) and (max_val_y > t[1]))
    ]


def get_diagonal_span(matrix, window, zero_tol):
    n = matrix.shape[0]
    lengths = [0] * n
    coords = [(0, 0)] * n

    for i in range(n):
        row = matrix[i]  # local view of row i
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
        coords[i] = (start * window, (end * window) + window)

    tuple_counts = Counter(coords)
    sorted_items = sorted(
        ((k, v) for k, v in tuple_counts.items() if v >= 3), key=lambda item: item[0][0]
    )
    return sorted_items


def merge_shared_boundaries(intervals, prefix, window, verbose=True):
    """
    intervals: list of ((start, end), count)
    returns: list of ((start+prefix, end+prefix), total_count)
    """
    out = []  # list of (x, y, count)

    chk_same = check_same_start_end  # expects: (pairs_list, val, window, is_start)
    chk_cont = check_new_contained  # expects: (pairs_list, x, y) checks if new sequence is smaller
    chk_span = (
        check_new_spans  # expects: (pairs_list, x, y) checks if new seq is bigger
    )

    for (x, y), count in intervals:
        # Base case, append to out if empty
        if len(out) == 0:
            out.append((x, y, count))
            continue
        else:
            # Initialize start and end window buffers
            start_range = (x - (window * 2), x + (window * 2))
            end_range = (y - (window * 2), y + (window * 2))
            process_entry(
                out=out,
                x=x,
                y=y,
                count=count,
                start_range=start_range,
                end_range=end_range,
                find_in_range_x=find_in_range_x,
                find_in_range_y=find_in_range_y,
                find_contained=find_contained,
                find_spanning=find_spanning,
                verbose=verbose,
            )
            print
    return [(x + prefix, y + prefix, cnt) for x, y, cnt in out]


def process_entry(
    out,
    x,
    y,
    count,
    start_range,
    end_range,
    find_in_range_x,
    find_in_range_y,
    find_contained,
    find_spanning,
    verbose=False,
):
    """
    Chooses the best prior entry to merge with (if any), then updates 'out'.
    Priority: match on both x & y (same index) > x > y > spanning > contained.
    Searches are done lazily to avoid unnecessary work.
    """

    idx = None
    reason = None

    # Try x match first (most common case) and only then check y for intersection
    c_x = find_in_range_x(data=out, min_val=start_range[0], max_val=start_range[1])
    if c_x:
        # Only compute y if x matched; try to find the same index for both
        c_y = find_in_range_y(data=out, min_val=end_range[0], max_val=end_range[1])
        if c_y:
            ix_x = {i for i, _ in c_x}
            ix_y = {i for i, _ in c_y}
            common = sorted(ix_x & ix_y)
            if common:
                idx = common[0]
                reason = "x & y"
            else:
                idx = c_x[0][0]
                reason = "x"
        else:
            idx = c_x[0][0]
            reason = "x"
    else:
        # No x match; try y
        c_y = find_in_range_y(data=out, min_val=end_range[0], max_val=end_range[1])
        if c_y:
            idx = c_y[0][0]
            reason = "y"
        else:
            # Only now attempt the more general/expensive checks
            c_s = find_spanning(
                data=out, min_val_x=x, max_val_x=x, min_val_y=y, max_val_y=y
            )
            if c_s:
                idx = c_s[0][0]
                reason = "spanning"
            else:
                c_c = find_contained(
                    data=out, min_val_x=x, max_val_x=x, min_val_y=y, max_val_y=y
                )
                if c_c:
                    idx = c_c[0][0]
                    reason = "contained"

    if idx is not None:
        if verbose:
            print(f"Matches previous entry on {reason}")
        _update_out(out, idx, x, y, count, verbose=verbose)
    else:
        out.append((x, y, count))


def sobel_spans(M, prefix):
    spans = []
    i = 0
    while i < len(M):
        # print(M[i])

        h = int(abs(M[i][3]))  # window length from height

        if h == 0 or h > 1900:
            i += 1
            continue

        end = min(i + h, len(M))  # clamp to array length
        if end <= i:  # safety, though h>0 makes this unlikely
            i += 1
            continue

        window = [abs(M[k][3]) for k in range(i, end)]
        med = median(window)

        if 0.75 <= med / h <= 1.25:
            # print("Yes")
            spans.append((i - 1, i + h + 1))
            i += h
        else:
            # print("No")
            i += 1
    spans = [(prefix * a, prefix * b) for (a, b) in spans]
    return spans


def sobel_with_diagonal_probes(M, thresh=0.7, min_thick=1):
    """
    Compute Sobel edges, display the binary edge map, and for each diagonal
    position (i, i) that is empty (False), draw a red vertical line that
    extends up and down until it hits a vertical run of True pixels with
    thickness >= min_thick.

    Parameters
    ----------
    M : 2D array
        Image/matrix to edge-detect.
    thresh : float
        Threshold on normalized Sobel magnitude to make binary_edges.
    min_thick : int
        Minimum contiguous thickness (in pixels) of a vertical edge to stop.
    figsize : tuple
        Matplotlib figure size.
    """
    # --- Sobel edges ---
    sobel_x = ndimage.sobel(M, axis=1)
    sobel_y = ndimage.sobel(M, axis=0)
    sobel_mag = np.hypot(sobel_x, sobel_y)
    max_val = np.max(sobel_mag)
    if max_val > 0:
        sobel_mag = sobel_mag / max_val
    binary_edges = sobel_mag > thresh

    H, W = binary_edges.shape
    N = min(H, W)

    def stop_y(y0, x, dy):
        """
        March from y0 in direction dy (+1 down, -1 up) until:
          - we hit image border, or
          - we encounter a vertical run of True pixels with length >= min_thick
            at column x, starting at the next step in the marching direction.
        Returns the last y BEFORE the blocking run/border.
        """
        y = y0
        while True:
            ny = y + dy
            if ny < 0 or ny >= H:
                return y  # hit border

            # Check if a vertical run with length >= min_thick begins at ny
            if dy > 0:
                end = min(ny + min_thick, H)
                if end - ny == min_thick and np.all(binary_edges[ny:end, x]):
                    return y
            else:  # dy < 0
                start = max(ny - (min_thick - 1), 0)
                if ny - start + 1 >= min_thick and np.all(
                    binary_edges[start : ny + 1, x]
                ):
                    return y

            y = ny  # keep marching

    # --- Plot base image and mask ---
    plt.figure(figsize=(8, 8))
    plt.imshow(binary_edges, cmap="gray_r", interpolation="nearest")
    plt.title("Sobel Edge Magnitude with Diagonal Probes")
    plt.colorbar(label="Edge (binary)")

    # --- For each diagonal position that is empty, drop a probe line ---
    for i in range(N):
        if not binary_edges[i, i]:  # "space on the diagonal"
            y_top = stop_y(i, i, dy=-1)
            y_bot = stop_y(i, i, dy=+1)
            # Draw the vertical line at x=i from y_top to y_bot
            plt.plot([i, i], [y_top, y_bot], "-", linewidth=1.5, color="red", alpha=0.9)

    plt.tight_layout()
    plt.show()
    ranges = []
    for i in range(N):
        if not binary_edges[i, i]:  # "space on the diagonal"
            y_top = stop_y(i, i, dy=-1)
            y_bot = stop_y(i, i, dy=+1)
            ranges.append((i, y_top, y_bot, y_top - y_bot))
    return ranges


def split_diagonal_attached(M):
    binary = M > 0
    vertical_structure = np.array([[0, 1, 0], [0, 1, 0], [0, 1, 0]], dtype=int)
    horizontal_structure = np.array([[0, 0, 0], [1, 1, 1], [0, 0, 0]], dtype=int)
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


def _update_out(out, index, x, y, count, verbose=False):
    old_x, old_y, old_count = out[index]
    min_x = min(old_x, x)
    max_y = max(old_y, y)
    # Keep old bounds if the new segment is "small" vs existing
    if count < old_count / 2:
        out[index] = (old_x, old_y, old_count + count)
    else:
        out[index] = (min_x, max_y, old_count + count)
    if verbose:
        print(f"Index: {index}")
        print(f"Updated to {out[index]}\n")
