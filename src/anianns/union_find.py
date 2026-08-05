from collections import defaultdict
import csv
from scipy import ndimage
import matplotlib.pyplot as plt
import numpy as np

try:
    import seaborn as sns
except ImportError:
    sns = None


# Function mapping colors to elements in the DSU
def assign_colors(items, palette_name="tab20"):
    """
    Assigns a unique color to each unique item in the input list.

    Parameters
    ----------
    items : list
        Input list (can have repeated elements)
    palette_name : str, optional
        Name of a seaborn/matplotlib palette (e.g. 'tab10', 'tab20', 'Set3', 'Paired', 'Spectral')

    Returns
    -------
    dict
        Mapping of unique items -> color (as RGB tuples)
    """
    unique_items = list(dict.fromkeys(items))  # preserves order & uniqueness
    n = len(unique_items)

    # Prefer seaborn palettes when available.
    if sns is not None:
        if n <= 20:
            palette = sns.color_palette(palette_name, n)
        else:
            # HUSL gives well-separated colors for large n.
            palette = sns.color_palette("husl", n)
    else:
        cmap_name = "tab20" if n <= 20 else "hsv"
        cmap = plt.get_cmap(cmap_name)
        palette = [cmap(i / max(1, n - 1))[:3] for i in range(n)]

    # Map each unique item to a color
    color_map = dict(zip(unique_items, palette))
    return color_map


def hashable_ntr_signature(value):
    """Convert composite NTRPrism metadata into an immutable color key."""
    if isinstance(value, np.ndarray):
        value = value.tolist()
    if isinstance(value, (list, tuple)):
        return tuple(hashable_ntr_signature(item) for item in value)
    if isinstance(value, dict):
        items = (
            (hashable_ntr_signature(key), hashable_ntr_signature(item))
            for key, item in value.items()
        )
        return tuple(sorted(items, key=repr))
    if isinstance(value, set):
        return frozenset(hashable_ntr_signature(item) for item in value)
    return value


def sobel_with_diagonal_probes2(M, thresh=0.6, min_thick=1, figsize=(8, 8)):
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
    plt.figure(figsize=figsize)
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
    return binary_edges


def sobel(M):
    sobel_x = ndimage.sobel(M, axis=1)
    sobel_y = ndimage.sobel(M, axis=0)
    sobel_mag = np.hypot(sobel_x, sobel_y)
    sobel_mag /= np.max(sobel_mag)  # Normalize
    binary_edges = sobel_mag > 0.6

    plt.imshow(binary_edges, cmap="gray_r")
    plt.title("Sobel Edge Magnitude")
    plt.colorbar()

    plt.show()
    return binary_edges


def sobel_edges(M, thresh=0.6, plot=True):
    """Return binary Sobel edge map."""
    sobel_x = ndimage.sobel(M, axis=1)
    sobel_y = ndimage.sobel(M, axis=0)
    sobel_mag = np.hypot(sobel_x, sobel_y)
    mx = np.max(sobel_mag)
    if mx > 0:
        sobel_mag = sobel_mag / mx
    binary_edges = sobel_mag > thresh
    if plot:
        dpi = 300
        figsize = (12, 10)
        plt.figure(figsize=figsize, dpi=dpi)
        plt.imshow(binary_edges, cmap="gray_r")
        plt.title("Sobel Edge Magnitude")
        plt.colorbar()

        plt.show()
    return binary_edges


def sobel_edges2(M, thresh=0.4, plot=True, highlight_ranges=None, offset=1.0):
    """Return binary Sobel edge map with optional highlighted regions."""
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import ndimage
    from matplotlib.patches import Rectangle

    sobel_x = ndimage.sobel(M, axis=1)
    sobel_y = ndimage.sobel(M, axis=0)
    sobel_mag = np.hypot(sobel_x, sobel_y)

    mx = np.max(sobel_mag)
    if mx > 0:
        sobel_mag = sobel_mag / mx

    binary_edges = sobel_mag > thresh

    if plot:
        dpi = 300
        figsize = (12, 10)

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        im = ax.imshow(binary_edges, cmap="gray_r")

        ax.set_title("Sobel Edge Magnitude")

        # --- optional: scale tick labels ---
        xticks = ax.get_xticks()
        yticks = ax.get_yticks()

        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
        ax.set_xticklabels([f"{x * offset:.2f}" for x in xticks])
        ax.set_yticklabels([f"{y * offset:.2f}" for y in yticks])

        # --- draw red highlight boxes ---
        if highlight_ranges is not None:
            for xmin, xmax, ymin, ymax in highlight_ranges:
                x0 = xmin / offset
                y0 = ymin / offset
                width = (xmax - xmin) / offset
                height = (ymax - ymin) / offset

                rect = Rectangle(
                    (x0, y0),
                    width,
                    height,
                    linewidth=2,
                    edgecolor="red",
                    facecolor="red",
                    alpha=0.3,
                )
                ax.add_patch(rect)

        plt.colorbar(im, ax=ax)
        plt.tight_layout()
        plt.show()

    return binary_edges


def sobel_edges3(M, thresh=0.4, plot=True, highlight_ranges=None, offset=1.0):
    """Return binary Sobel edge map with optional highlighted regions."""
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import ndimage
    from matplotlib.patches import Rectangle

    sobel_x = ndimage.sobel(M, axis=1)
    sobel_y = ndimage.sobel(M, axis=0)
    sobel_mag = np.hypot(sobel_x, sobel_y)

    mx = np.max(sobel_mag)
    if mx > 0:
        sobel_mag = sobel_mag / mx

    binary_edges = sobel_mag > thresh

    if plot:
        dpi = 300
        figsize = (12, 10)

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        im = ax.imshow(binary_edges, cmap="gray_r")

        ax.set_title("Sobel Edge Magnitude")

        # --- optional: scale tick labels ---
        xticks = ax.get_xticks()
        yticks = ax.get_yticks()

        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
        ax.set_xticklabels([f"{x * offset:.2f}" for x in xticks])
        ax.set_yticklabels([f"{y * offset:.2f}" for y in yticks])

        # --- draw red highlight boxes ---
        if highlight_ranges is not None:
            for xmin, xmax, ymin, ymax in highlight_ranges:
                x0 = xmin / offset
                y0 = ymin / offset
                width = (xmax - xmin) / offset
                height = (ymax - ymin) / offset

                rect = Rectangle(
                    (x0, y0),
                    width,
                    height,
                    linewidth=2,
                    edgecolor="green",
                    facecolor="green",
                    alpha=0.3,
                )
                ax.add_patch(rect)

        plt.colorbar(im, ax=ax)
        plt.tight_layout()
        plt.show()

    return binary_edges


def _march_until_hit(binary_edges, y0, x0, dy, dx):
    """
    March from (y0, x0) in direction (dy, dx) until the first True pixel.

    Returns
    -------
    tuple
        ((y, x), hit_found) where (y, x) is either the first hit coordinate
        or the last in-bounds coordinate reached if no hit is found.
    """
    h, w = binary_edges.shape
    y, x = y0, x0

    while True:
        ny = y + dy
        nx = x + dx
        if ny < 0 or ny >= h or nx < 0 or nx >= w:
            return (y, x), False
        if binary_edges[ny, nx]:
            return (ny, nx), True
        y, x = ny, nx


def probe_diagonal_rays(binary_edges):
    """
    For each diagonal position (i, i), shoot rays up/down/left/right and
    report the first hit in each direction.

    If no hit is found in a direction, the returned coordinate is the last
    in-bounds position reached along that ray.
    """
    h, w = binary_edges.shape
    n = min(h, w)
    reports = []

    for i in range(n):
        up, up_hit = _march_until_hit(binary_edges, i, i, -1, 0)
        down, down_hit = _march_until_hit(binary_edges, i, i, 1, 0)
        left, left_hit = _march_until_hit(binary_edges, i, i, 0, -1)
        right, right_hit = _march_until_hit(binary_edges, i, i, 0, 1)

        reports.append(
            {
                "diag_index": i,
                "origin": (i, i),
                "up": up,
                "down": down,
                "left": left,
                "right": right,
                "up_hit": up_hit,
                "down_hit": down_hit,
                "left_hit": left_hit,
                "right_hit": right_hit,
            }
        )

    return reports


def find_diagonal_bisecting_squares(
    M, thresh=0.6, min_span=1, square_tol=0, plot=False, edge_map=None
):
    """
    Find square-like boxes that bisect the diagonal by ray-casting from each
    diagonal coordinate.

    For each diagonal point, the candidate box is defined by the first hit
    upward, downward, leftward, and rightward. A box is reported only when all
    four rays hit an edge and the resulting height/width differ by at most
    ``square_tol``.
    """
    if edge_map is not None:
        edges = np.asarray(edge_map, dtype=bool)
    else:
        arr = np.asarray(M)
        if arr.dtype == np.bool_:
            edges = arr
        else:
            edges = sobel_edges(arr, thresh=thresh, plot=False)

    reports = probe_diagonal_rays(edges)
    squares = []
    seen = set()

    for report in reports:
        if not all(
            report[key] for key in ("up_hit", "down_hit", "left_hit", "right_hit")
        ):
            continue

        diag_index = report["diag_index"]
        top = report["up"][0]
        bottom = report["down"][0]
        left = report["left"][1]
        right = report["right"][1]

        height = bottom - top
        width = right - left

        if height < min_span or width < min_span:
            continue
        if abs(height - width) > square_tol:
            continue
        if not (top < diag_index < bottom and left < diag_index < right):
            continue

        square = (top, left, bottom, right)
        if square not in seen:
            seen.add(square)
            squares.append(square)

    if plot:
        plt.figure(figsize=(8, 8))
        plt.imshow(edges, cmap="gray_r", interpolation="nearest")
        ax = plt.gca()

        for top, left, bottom, right in squares:
            ax.add_patch(
                plt.Rectangle(
                    (left, top),
                    right - left,
                    bottom - top,
                    edgecolor="cyan",
                    facecolor="none",
                    lw=2,
                )
            )

        plt.title("Diagonal-Bisecting Squares")
        plt.tight_layout()
        plt.show()

    return reports, squares


def find_offdiag_rectangles(
    M, thresh=0.4, band=5, connectivity=2, min_height=1, min_width=1, plot=True
):
    """
    Find bounding boxes for connected components that do NOT intersect
    a diagonal band |row-col| <= band. Returns list of (y0, x0, y1, x1).
    """
    edges = sobel_edges(M, thresh=thresh, plot=plot)
    H, W = edges.shape

    # Connected components on the full edge map
    structure = ndimage.generate_binary_structure(
        2, connectivity
    )  # 4-conn if 1, 8-conn if 2
    labeled, num = ndimage.label(edges, structure=structure)
    slices = ndimage.find_objects(labeled)

    rects = []
    if slices is not None:
        for label_id, sl in enumerate(slices, start=1):
            if sl is None:
                continue
            sy, sx = sl
            y0, y1 = sy.start, sy.stop
            x0, x1 = sx.start, sx.stop

            # Optional size filter
            if (y1 - y0) < min_height or (x1 - x0) < min_width:
                continue

            # Gather pixels in this component
            comp_mask = labeled[sy, sx] == label_id
            if not np.any(comp_mask):
                continue
            # Compute deltas (row - col) in the subwindow
            rr, cc = np.indices(comp_mask.shape)
            rr += y0
            cc += x0
            deltas = (rr - cc)[comp_mask]

            # "Bisect" criterion: crosses both sides of diagonal band
            # i.e., has pixels with delta <= -band and delta >= +band
            if deltas.min() <= -band and deltas.max() >= +band:
                rects.append((y0, x0, y1, x1))

    if plot:
        plt.figure(figsize=(8, 8))
        plt.imshow(edges, cmap="gray_r", interpolation="nearest")

        # Visualize the diagonal band to show where "bisecting" occurs
        rr, cc = np.indices((H, W))
        band_mask = np.abs(rr - cc) <= band
        band_vis = np.full_like(edges, np.nan, dtype=float)
        band_vis[band_mask] = 1.0
        plt.imshow(band_vis, cmap="Reds", alpha=0.25)  # diagonal band overlay

        # Draw rectangles for bisecting components
        ax = plt.gca()
        for yy0, xx0, yy1, xx1 in rects:
            ax.add_patch(
                plt.Rectangle(
                    (xx0, yy0),
                    xx1 - xx0,
                    yy1 - yy0,
                    edgecolor="lime",
                    facecolor="none",
                    lw=2,
                )
            )
        plt.title("Diagonal-bisecting component rectangles")
        plt.tight_layout()
        plt.show()

    return rects


class DSU:
    def __init__(self, n=0):
        self.parent = list(range(n))
        self.size = [1] * n

    def _grow_to(self, n):
        # Ensure DSU has capacity for n items total
        while len(self.parent) < n:
            i = len(self.parent)
            self.parent.append(i)
            self.size.append(1)

    def find(self, x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        if self.size[ra] < self.size[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        self.size[ra] += self.size[rb]
        return True


class SatelliteDSU:
    """Disjoint-set collection of refined satellites linked by distal hits."""

    def __init__(self):
        self.satellites = []
        self.dsu = DSU(0)
        self._key_to_index = {}
        self._neighbors = defaultdict(set)

    def add_satellite(self, chrom, start, end, **metadata):
        start = int(start)
        end = int(end)
        if end <= start:
            raise ValueError("satellite end must be greater than start")
        key = (str(chrom), start, end)
        existing = self._key_to_index.get(key)
        if existing is not None:
            for name, value in metadata.items():
                if value is not None:
                    self.satellites[existing][name] = value
            return existing

        index = len(self.satellites)
        satellite = {"chrom": key[0], "start": start, "end": end}
        satellite.update(metadata)
        self.satellites.append(satellite)
        self._key_to_index[key] = index
        self.dsu._grow_to(index + 1)
        return index

    @staticmethod
    def _reciprocal_overlap(start1, end1, start2, end2):
        overlap = max(0, min(end1, end2) - max(start1, start2))
        if overlap == 0:
            return 0.0
        return min(overlap / (end1 - start1), overlap / (end2 - start2))

    def find_satellite(self, chrom, start, end, minimum_overlap=0.5):
        """Resolve a distal endpoint to its best refined satellite node."""
        chrom = str(chrom)
        start = int(start)
        end = int(end)
        exact = self._key_to_index.get((chrom, start, end))
        if exact is not None:
            return exact
        best_index = None
        best_overlap = 0.0
        for index, satellite in enumerate(self.satellites):
            if satellite["chrom"] != chrom:
                continue
            overlap = self._reciprocal_overlap(
                start,
                end,
                satellite["start"],
                satellite["end"],
            )
            if overlap > best_overlap:
                best_index = index
                best_overlap = overlap
        return best_index if best_overlap >= minimum_overlap else None

    def union(self, first_index, second_index):
        if first_index == second_index:
            return False
        self._neighbors[first_index].add(second_index)
        self._neighbors[second_index].add(first_index)
        return self.dsu.union(first_index, second_index)

    def union_by_coordinates(
        self,
        chrom1,
        start1,
        end1,
        chrom2,
        start2,
        end2,
        minimum_overlap=0.5,
    ):
        """Union the two refined satellites matching a distal relationship."""
        first = self.find_satellite(chrom1, start1, end1, minimum_overlap)
        second = self.find_satellite(chrom2, start2, end2, minimum_overlap)
        if first is None or second is None:
            return False
        self.union(first, second)
        return True

    def component_rows(self):
        """Return deterministic rows describing every DSU member."""
        buckets = defaultdict(list)
        for index in range(len(self.satellites)):
            buckets[self.dsu.find(index)].append(index)

        groups = sorted(
            buckets.values(),
            key=lambda members: min(
                (
                    self.satellites[index]["chrom"],
                    self.satellites[index]["start"],
                    self.satellites[index]["end"],
                )
                for index in members
            ),
        )
        component_for_index = {}
        component_size = {}
        for component_number, members in enumerate(groups, start=1):
            component_id = f"satellite_component_{component_number:04d}"
            for index in members:
                component_for_index[index] = component_id
                component_size[index] = len(members)

        sorted_indices = sorted(
            range(len(self.satellites)),
            key=lambda index: (
                self.satellites[index]["chrom"],
                self.satellites[index]["start"],
                self.satellites[index]["end"],
            ),
        )
        rows = []
        for satellite_number, index in enumerate(sorted_indices, start=1):
            satellite = self.satellites[index]
            rows.append(
                {
                    "satellite_id": f"satellite_{satellite_number:06d}",
                    "component_id": component_for_index[index],
                    "component_size": component_size[index],
                    "chrom": satellite["chrom"],
                    "start": satellite["start"],
                    "end": satellite["end"],
                    "name": satellite.get("name"),
                    "monomer": satellite.get("monomer"),
                    "periodicity": satellite.get("periodicity"),
                    "is_hor": satellite.get("is_hor"),
                    "direct_link_count": len(self._neighbors[index]),
                }
            )
        return rows

    def item_rgb_for_annotations(
        self,
        chrom,
        starts,
        ends,
        monomers,
        periodicities,
        hor_flags,
    ):
        """Color annotations by DSU component and exact NTRPrism signature."""
        values = list(zip(starts, ends, monomers, periodicities, hor_flags))
        color_keys = []
        for row_number, (start, end, monomer, periodicity, is_hor) in enumerate(values):
            index = self.find_satellite(chrom, start, end)
            if index is None or monomer in (None, 0):
                # Missing NTRPrism evidence must not cause unrelated calls to
                # appear matched simply because both values are absent.
                color_keys.append(("unmatched", str(chrom), row_number))
                continue
            color_keys.append(
                (
                    "matched",
                    self.dsu.find(index),
                    hashable_ntr_signature(monomer),
                    hashable_ntr_signature(periodicity),
                    bool(is_hor),
                )
            )

        color_map = assign_colors(color_keys)
        return [
            ",".join(
                str(min(255, max(0, int(round(channel * 255)))))
                for channel in color_map[key]
            )
            for key in color_keys
        ]

    def write_tsv(self, output_path):
        """Serialize all nodes and their final connected-component IDs."""
        fieldnames = (
            "satellite_id",
            "component_id",
            "component_size",
            "chrom",
            "start",
            "end",
            "name",
            "monomer",
            "periodicity",
            "is_hor",
            "direct_link_count",
        )
        with open(output_path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(self.component_rows())

    def write_text(self, output_path):
        """Write a human-readable summary grouped by DSU component."""
        rows = self.component_rows()
        components = defaultdict(list)
        for row in rows:
            components[row["component_id"]].append(row)
        linked_component_count = sum(
            len(members) > 1 for members in components.values()
        )

        with open(output_path, "w") as handle:
            handle.write("Satellite DSU Results\n")
            handle.write("=====================\n")
            handle.write(f"Total satellites: {len(rows)}\n")
            handle.write(f"Total components: {len(components)}\n")
            handle.write(f"Linked components: {linked_component_count}\n\n")

            for component_id, members in components.items():
                handle.write(f"[{component_id}] size={len(members)}\n")
                for row in members:
                    name = row["name"] if row["name"] is not None else "Unknown"
                    monomer = (
                        row["monomer"] if row["monomer"] is not None else "Unknown"
                    )
                    handle.write(
                        f"  {row['satellite_id']} "
                        f"{row['chrom']}:{row['start']}-{row['end']} "
                        f"name={name} monomer={monomer} "
                        f"direct_links={row['direct_link_count']}\n"
                    )
                handle.write("\n")


class TupleDSU:
    """
    Stores:
      - 'point' items: (name, x, y, color)
      - 'edge' items: (x0, y0, x1, y1)
    Unions happen when coordinates match:
      - point matches point at (x,y)
      - edge matches anything touching (x0,y0) and/or (x1,y1)
    """

    def __init__(self):
        self.items = []  # heterogeneous: dicts with 'kind' field
        self.dsu = DSU(0)
        self.xy_map = defaultdict(
            list
        )  # (x,y) -> list of indices registered at that coordinate

    def _add_item(self, item, coords_to_register):
        """
        Internal helper:
          - append item
          - grow DSU
          - register item under each coord in coords_to_register
          - union with all previously registered indices at those coords
        """
        idx = len(self.items)
        self.items.append(item)
        self.dsu._grow_to(idx + 1)

        for coord in coords_to_register:
            # union with everything already at this coordinate
            for j in self.xy_map[coord]:
                self.dsu.union(idx, j)
            # then register this new index for future matches
            self.xy_map[coord].append(idx)

        return idx

    # --- Public APIs ---

    def add_point(self, name, x, y, color):
        """
        Add a point-tuple (name, x, y, color).
        Unions with any items already registered at (x,y).
        """
        item = {"kind": "point", "name": name, "x": x, "y": y, "color": color}
        return self._add_item(item, coords_to_register=[(x, y)])

    def add_edge(self, x0, y0, x1, y1):
        """
        Add an edge (x0, y0, x1, y1).
        Unions with anything at (x0, y0) and anything at (x1, y1).
        If both endpoints match existing items, this edge will connect those groups.
        """
        item = {"kind": "edge", "x0": x0, "y0": y0, "x1": x1, "y1": y1}
        return self._add_item(item, coords_to_register=[(x0, y0), (x1, y1)])

    def groups(self):
        """
        Returns groups as lists of stored items (points and/or edges).
        """
        from collections import defaultdict

        buckets = defaultdict(list)
        for i, it in enumerate(self.items):
            root = self.dsu.find(i)
            buckets[root].append(it)
        return list(buckets.values())
