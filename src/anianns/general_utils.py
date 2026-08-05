import os
import json
from typing import Union, List, Dict
from anianns.const import BED_COLUMNS

# from anianns.kmer_utils import generateKmersFromFasta, generateKmersFromFastaForwardOnly
from itertools import islice
import polars as pl
import pysam
from collections import Counter
import numpy as np
import re
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import csv
from io import StringIO
from urllib.parse import quote

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans"]


def clean_genomic_ticks(start, end, target_intervals=4):
    """Return readable genomic tick values while retaining exact band edges."""
    locator = MaxNLocator(
        nbins=target_intervals,
        steps=[1, 2, 2.5, 5, 10],
    )
    ticks = locator.tick_values(start, end)
    tolerance = max(abs(end - start), 1.0) * 1e-9
    ticks = ticks[(ticks >= start - tolerance) & (ticks <= end + tolerance)]
    if ticks.size == 0 or not np.isclose(ticks[0], start):
        ticks = np.insert(ticks, 0, start)
    else:
        ticks[0] = start
    if not np.isclose(ticks[-1], end):
        ticks = np.append(ticks, end)
    else:
        ticks[-1] = end
    return ticks

from anianns.kmer_utils import (
    generate_kmers_from_fasta,
    generate_kmers_from_fasta_forward_only,
)


def add_prefix_to_tuples(data, band_height, w):
    prefix_amount = band_height * (w - 1)
    return [(start + prefix_amount, end + prefix_amount) for start, end in data]


def calculate_distances(numbers):
    indices_map = {}
    distances = []

    # Iterate through the list and store indices of each number
    for idx, num in enumerate(numbers):
        if num in indices_map:
            # Calculate the distance betwe en current and previous occurrence
            prev_idx = indices_map[num][-1]
            distance = idx - prev_idx
            distances.append(distance)
            # Update the list of indices for this number
            indices_map[num].append(idx)
        else:
            # If the number is seen for the first time, just store its index
            indices_map[num] = [idx]

    return distances


def check_bed_vs_indexed_fasta(
    bed_dfs: List[pl.DataFrame], fasta_paths: Union[str, List[str]]
) -> None:
    fasta_chroms = get_fasta_indexed_chroms(fasta_paths)

    for i, df in enumerate(bed_dfs):
        if "chrom" not in df.columns:
            print(f"[WARNING] BED file {i} missing 'chrom' column — skipping.")
            continue

        bed_chroms = set(df.select("chrom").unique().to_series(0).to_list())
        missing = bed_chroms - fasta_chroms
        if missing:
            print(
                f"Input fastas missing: {sorted(missing)}. Check 'chr' spelling or contents in .fai file...\n"
            )
            return False
    return True


def convert_dataframe_format(df: pl.DataFrame, format: str) -> str:
    """
    Convert a Polars DataFrame with BED-like schema into a given format.

    Parameters
    ----------
    df : pl.DataFrame
        The input dataframe with columns:
        chrom, start, end, name, score, strand, thickStart, thickEnd, itemRgb
    format : str
        Output format: one of ["gtf", "gff", "csv", "tsv", "json"]

    Returns
    -------
    str
        The converted data as a string in the requested format.
    """
    format = format.lower()

    chrom_column = "chrom" if "chrom" in df.columns else "#chrom"
    required = {chrom_column, "start", "end", "name", "score", "strand"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(
            "Input annotation dataframe is missing required column(s): "
            + ", ".join(sorted(missing))
        )

    def gtf_escape(value):
        return (
            str(value)
            .replace("\\", "\\\\")
            .replace('"', '\\"')
            .replace("\t", " ")
            .replace("\r", " ")
            .replace("\n", " ")
        )

    # --- BED → GTF conversion ---
    if format == "gtf":
        output = StringIO()
        for index, row in enumerate(df.iter_rows(named=True), start=1):
            start, end = int(row["start"]), int(row["end"])
            if end <= start:
                continue
            repeat_id = f"anianns_{index:06d}"
            repeat_name = gtf_escape(row["name"] or "Unclassified Repeat")
            monomer = "." if row["score"] is None else str(row["score"])
            attributes = (
                f'gene_id "{repeat_id}"; '
                f'transcript_id "{repeat_id}"; '
                f'repeat_name "{repeat_name}"; '
                f'monomer_length "{gtf_escape(monomer)}";'
            )
            fields = (
                row[chrom_column],
                "AniAnns",
                "tandem_repeat",
                start + 1,
                end,
                ".",
                row["strand"] or ".",
                ".",
                attributes,
            )
            output.write("\t".join(map(str, fields)) + "\n")
        return output.getvalue()

    # --- BED → GFF conversion ---
    elif format == "gff":
        output = StringIO()
        output.write("##gff-version 3\n")
        for index, row in enumerate(df.iter_rows(named=True), start=1):
            start, end = int(row["start"]), int(row["end"])
            if end <= start:
                continue
            repeat_id = f"anianns_{index:06d}"
            repeat_name = quote(str(row["name"] or "Unclassified Repeat"), safe="")
            monomer = (
                "."
                if row["score"] is None
                else quote(str(row["score"]), safe="")
            )
            attributes = (
                f"ID={repeat_id};Name={repeat_name};monomer_length={monomer}"
            )
            fields = (
                row[chrom_column],
                "AniAnns",
                "tandem_repeat",
                start + 1,
                end,
                ".",
                row["strand"] or ".",
                ".",
                attributes,
            )
            output.write("\t".join(map(str, fields)) + "\n")
        return output.getvalue()

    # --- CSV ---
    elif format == "csv":
        return df.write_csv(None)

    # --- TSV ---
    elif format == "tsv":
        return df.write_csv(None, separator="\t")

    # --- JSON ---
    elif format == "json":
        return df.write_json()

    else:
        raise ValueError(
            "Invalid format. Choose from: 'gtf', 'gff', 'csv', 'tsv', 'json'."
        )


def define_bounds(seq_name):
    """Extract chromosome and region from seq_name.

    Supports:
    - Standard format: "chrY:50-3000"
    - Extended format: "HG002_chr13_MATERNAL:1-4000000:1000000-3000000"
      (keeps the last range).
    """
    # Match chromosome + one or more ranges separated by colons
    region_pattern = r"^([a-zA-Z0-9_]+)(?::(\d+-\d+))+"
    match = re.match(region_pattern, seq_name)
    if match:
        # Extract chromosome, lower bound, and upper bound
        chrom = match.group(1)
        ranges = re.findall(r"(\d+)-(\d+)", seq_name)
        if ranges:
            # Take the last range
            lower_bound, upper_bound = map(int, ranges[-1])
            return chrom, lower_bound, upper_bound

    # No match
    return None


def extract_region(fasta_file, chr, region_start, region_end):
    """
    Extract a sequence from a FASTA file using 1-based coordinates (inclusive).

    Args:
        fasta_file: Path to an indexed FASTA, or an open ``pysam.FastaFile``.
        chr (str): Chromosome or contig name.
        region_start (int): 1-based start coordinate.
        region_end (int): 1-based end coordinate (inclusive).

    Returns:
        str or None: Extracted DNA sequence, or None if an error occurred.
    """
    fasta = None
    owns_handle = False
    try:
        if region_start < 1:
            region_start = 1
        if hasattr(fasta_file, "fetch"):
            fasta = fasta_file
        else:
            fasta = pysam.FastaFile(fasta_file)
            owns_handle = True
        sequence = fasta.fetch(chr, region_start, region_end)
        return sequence
    except Exception as e:
        print(
            f"Error fetching region {chr}:{region_start}-{region_end} from {fasta_file}\n"
            f"Details: {e}\n"
        )
        return None
    finally:
        if owns_handle and fasta is not None:
            fasta.close()


def extract_regions_by_name(
    df: pl.DataFrame,
    fasta_files: Union[str, List[str]],
    k: int,
    verbose: bool,
) -> Dict[str, set]:
    """
    Given a DataFrame with columns 'chrom', 'start', 'end', and 'name',
    and one or more FASTA file paths (string or list of strings),
    extract each region from all provided FASTA files, generate k-mers,
    and return a dict mapping each lowercase name to its combined k-mer set,
    printing the FASTA file each sequence is retrieved from.
    """
    # Ensure fasta_files is a list
    if isinstance(fasta_files, str):
        fasta_files = [fasta_files]

    kmer_dict: Dict[str, set] = {}
    unique_names = df["name"].unique().to_list()

    for name in unique_names:
        sub_df = df.filter(pl.col("name") == name)
        print(f"Retrieving k-mers for '{name}'...")

        kmer_set = set()
        for row in sub_df.iter_rows(named=True):
            chrom, start, end = row["chrom"], row["start"], row["end"]
            # Attempt to extract from each FASTA until found
            for fasta in fasta_files:
                seq = extract_region(fasta, chrom, start, end)
                if seq:
                    if verbose:
                        print(
                            f"  Retrieved region {chrom}:{start}-{end} from: '{fasta}'"
                        )
                    kmers = generate_kmers_from_fasta(seq, k, quiet=True)
                    kmer_set.update(kmers)
                    break
            else:
                print(
                    f"  Warning: region {chrom}:{start}-{end} not found in any provided FASTA."
                )

        kmer_dict[name.lower()] = kmer_set

    return kmer_dict


def extract_histograms_by_name(
    df: pl.DataFrame, fasta_files: Union[str, List[str]], k: int, verbose: bool
) -> Dict[str, List]:
    """
    Given a Polars DataFrame with columns 'chrom', 'start', 'end', and 'name',
    and one or more FASTA file paths (string or list of strings),
    extract each region, generate forward-only k-mers,
    compute distance histograms, and return a dict mapping each lowercase
    name to its list of top-10 frequent distance histograms,
    printing which FASTA each sequence came from.
    """
    # Normalize FASTA input to a list
    if isinstance(fasta_files, str):
        fasta_files = [fasta_files]

    cdf_dict: Dict[str, List] = {}
    unique_names = df["name"].unique().to_list()

    for name in unique_names:
        sub_df = df.filter(pl.col("name") == name)
        print(f"Retrieving k-mer histograms for '{name}'...")

        kmer_histograms: List = []
        for row in sub_df.iter_rows(named=True):
            chrom, start, end = row["chrom"], row["start"], row["end"]
            # Try each FASTA file until a sequence is found
            for fasta in fasta_files:
                seq = extract_region(fasta, chrom, start, end)
                if seq:
                    if verbose:
                        print(
                            f"  Retrieved region {chrom}:{start}-{end} from: '{fasta}'"
                        )
                    size = end - start
                    kmers = generate_kmers_from_fasta_forward_only(seq, k, quiet=True)
                    kmer_list = list(islice(kmers, 1, size))
                    histo_list = calculate_distances(kmer_list)
                    top_dists = top_n_frequent_distances(histo_list, 10)
                    if verbose:
                        print(name.lower(), top_dists)
                    kmer_histograms.append(top_dists)
                    break
            else:
                print(
                    f"  Warning: region {chrom}:{start}-{end} not found in any provided FASTA."
                )

        cdf_dict[name.lower()] = kmer_histograms

    return cdf_dict


def get_fasta_indexed_chroms(fasta_paths: Union[str, List[str]]) -> set:
    if isinstance(fasta_paths, str):
        fasta_paths = [fasta_paths]

    chroms = set()
    for path in fasta_paths:
        try:
            with pysam.FastaFile(path) as fasta:
                chroms.update(fasta.references)
        except Exception as e:
            print(f"[ERROR] Failed to open or index FASTA file: {path}\n{e}")
    return chroms


def get_input_headers(filename: List) -> List:
    header_list = []
    for file in filename:
        try:
            seq_list = []
            seq = pysam.FastaFile(file)
            for seq_id in seq.references:
                seq_list.append(seq_id)
            header_list.append((file, seq_list))
        except OSError:
            seq = None

    return header_list


def merge_close_values(pairs, tolerance=1):
    # sort by value first for correct merging
    pairs = sorted(pairs, key=lambda x: x[0])

    merged = []

    for value, count in pairs:
        if not merged:
            merged.append(
                {"values": [(value, count)], "total_count": count, "rep_value": value}
            )
            continue

        last = merged[-1]

        if abs(last["rep_value"] - value) <= tolerance:
            last["values"].append((value, count))
            last["total_count"] += count
            last["rep_value"] = max(last["values"], key=lambda x: x[1])[0]
        else:
            merged.append(
                {"values": [(value, count)], "total_count": count, "rep_value": value}
            )

    # 🔑 sort final output by total_count (descending)
    result = [(g["rep_value"], g["total_count"]) for g in merged]
    result.sort(key=lambda x: x[1], reverse=True)

    return result


def plot_matrix(
    matrix,
    title="Matrix Plot",
    cmap="gray_r",
    show_colorbar=True,
    dpi=72,
    figsize=(6, 5),
    save_path=None,
    offset=1.0,
    highlight_ranges=None,  # list of (xmin, xmax, ymin, ymax)
    diagonal_ranges=None,  # list of (start, end)
    edge_overlay=None,
    edge_only=False,
    colorbar_label="Value",
    colorbar_pad=0.04,
    reserve_colorbar_space=False,
    coordinate_origin=None,
    coordinate_end=None,
    coordinate_units="bp",
    white_below=None,
    vmin=86,
    vmax=100,
    aspect="auto",
    legend_outside=False,
):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    if not isinstance(matrix, np.ndarray):
        raise TypeError("Input must be a NumPy array")
    if matrix.ndim != 2:
        raise ValueError("Input must be a 2D matrix")
    if edge_overlay is not None:
        edge_overlay = np.asarray(edge_overlay, dtype=np.float32)
        if edge_overlay.shape != matrix.shape:
            raise ValueError("Sobel edge overlay must match the matrix shape")
        if not np.all(np.isfinite(edge_overlay)):
            raise ValueError("Sobel edge overlay must contain finite values")
        edge_overlay = np.clip(edge_overlay, 0.0, 1.0)
    if edge_only and edge_overlay is None:
        raise ValueError("edge-only plotting requires a Sobel edge overlay")

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    if cmap == "spectral_11":
        cmap = plt.get_cmap("Spectral", 11)
    elif cmap == "spectral_11_r":
        cmap = plt.get_cmap("Spectral_r", 11)

    displayed_matrix = np.zeros_like(matrix) if edge_only else matrix
    if white_below is not None:
        displayed_matrix = np.ma.masked_less(displayed_matrix, white_below)
        cmap = cmap.copy() if hasattr(cmap, "copy") else cmap
        cmap.set_bad("white")
    im = ax.imshow(
        displayed_matrix,
        cmap=cmap,
        aspect=aspect,
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
    )
    if edge_overlay is not None and np.any(edge_overlay):
        from matplotlib.colors import LinearSegmentedColormap

        visible_edges = np.ma.masked_where(edge_overlay <= 0.05, edge_overlay)
        edge_cmap = LinearSegmentedColormap.from_list(
            "sobel_cyan",
            [(0.0, 0.75, 1.0, 0.0), (0.0, 0.75, 1.0, 1.0)],
        )
        ax.imshow(
            visible_edges,
            cmap=edge_cmap,
            aspect=aspect,
            interpolation="nearest",
            vmin=0,
            vmax=1,
            zorder=2,
        )
    if aspect == "equal":
        ax.set_box_aspect(1)

    ax.set_title(title)
    if coordinate_origin is None:
        ax.set_xlabel("Columns")
        ax.set_ylabel("Rows")
    else:
        axis_label = f"Genomic position ({coordinate_units})"
        ax.set_xlabel(axis_label)
        ax.set_ylabel(axis_label)

    # --- scale tick labels only ---
    if coordinate_origin is not None:
        divisor = 1_000_000 if coordinate_units == "Mbp" else 1
        x_end = coordinate_end or coordinate_origin + matrix.shape[1] * offset
        y_end = coordinate_end or coordinate_origin + matrix.shape[0] * offset
        x_coordinates = clean_genomic_ticks(
            coordinate_origin / divisor, x_end / divisor
        )
        y_coordinates = clean_genomic_ticks(
            coordinate_origin / divisor, y_end / divisor
        )
        x_fraction = (x_coordinates * divisor - coordinate_origin) / (
            x_end - coordinate_origin
        )
        y_fraction = (y_coordinates * divisor - coordinate_origin) / (
            y_end - coordinate_origin
        )
        xticks = -0.5 + (x_fraction * matrix.shape[1])
        yticks = -0.5 + (y_fraction * matrix.shape[0])
    elif aspect == "equal":
        xticks = np.linspace(-0.5, max(-0.5, matrix.shape[1] - 0.5), 5)
        yticks = np.linspace(-0.5, max(-0.5, matrix.shape[0] - 0.5), 5)
    else:
        xticks = ax.get_xticks()
        yticks = ax.get_yticks()

    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    if coordinate_origin is not None:
        ax.set_xticklabels([f"{value:g}" for value in x_coordinates])
        ax.set_yticklabels([f"{value:g}" for value in y_coordinates])
    elif aspect == "equal":
        ax.set_xticklabels([f"{x * offset:,.0f}" for x in xticks])
        ax.set_yticklabels([f"{y * offset:,.0f}" for y in yticks])
    else:
        ax.set_xticklabels([f"{x * offset:.2f}" for x in xticks])
        ax.set_yticklabels([f"{y * offset:.2f}" for y in yticks])

    # --- detected diagonal satellites ---
    if diagonal_ranges is not None:
        for start, end in diagonal_ranges:
            start_index = start / offset
            end_index = end / offset
            rect = Rectangle(
                (start_index, start_index),
                end_index - start_index,
                end_index - start_index,
                linewidth=2,
                edgecolor="#00A651",
                facecolor="none",
                zorder=3,
            )
            ax.add_patch(rect)

    # --- distal highlight regions ---
    if highlight_ranges is not None:
        for xmin, xmax, ymin, ymax in highlight_ranges:
            # convert scaled values → indices
            x0 = xmin / offset
            x1 = xmax / offset
            y0 = ymin / offset
            y1 = ymax / offset

            rect = Rectangle(
                (x0, y0),
                x1 - x0,
                y1 - y0,
                linewidth=2,
                edgecolor="red",
                facecolor="red",
                alpha=0.3,
                zorder=3,
            )
            ax.add_patch(rect)

    legend_handles = []
    if edge_overlay is not None and np.any(edge_overlay):
        from matplotlib.lines import Line2D

        legend_handles.append(
            Line2D([0], [0], color="#00BFFF", linewidth=2, label="Sobel edges")
        )
    if diagonal_ranges:
        legend_handles.append(
            Rectangle(
                (0, 0),
                1,
                1,
                edgecolor="#00A651",
                facecolor="none",
                label="Detected satellite",
            )
        )
    if highlight_ranges:
        legend_handles.append(
            Rectangle(
                (0, 0),
                1,
                1,
                edgecolor="red",
                facecolor="red",
                alpha=0.3,
                label="Distal link",
            )
        )
    if legend_handles:
        if legend_outside:
            fig.legend(
                handles=legend_handles,
                loc="center left",
                bbox_to_anchor=(0.78, 0.5),
                borderaxespad=0,
                fontsize=6,
            )
        else:
            ax.legend(handles=legend_handles, loc="lower right", fontsize=6)

        """# --- cross-region (off-diagonal) rectangles ---
    if highlight_ranges is not None:
        for i in range(len(highlight_ranges)):
            for j in range(i + 1, len(highlight_ranges)):  # avoid duplicates

                x1_min, x1_max, _, _ = highlight_ranges[i]
                x2_min, x2_max, _, _ = highlight_ranges[j]

                a1, b1 = x1_min, x1_max
                a2, b2 = x2_min, x2_max

                # ✅ ONLY allow strictly separated ranges
                if not (b1 < a2 - 5 or b2 < a1 - 5):
                    continue  # skip overlaps AND touching

                # --- rectangle: (j on x-axis, i on y-axis)
                x0 = a2 / offset
                y0 = a1 / offset
                width = (b2 - a2) / offset
                height = (b1 - a1) / offset

                rect = Rectangle(
                    (x0, y0),
                    width,
                    height,
                    linewidth=1.5,
                    edgecolor="green",
                    facecolor="green",
                    alpha=0.2,
                )
                ax.add_patch(rect)

                # --- optional symmetric rectangle
                x0_sym = a1 / offset
                y0_sym = a2 / offset

                rect_sym = Rectangle(
                    (x0_sym, y0_sym),
                    (b1 - a1) / offset,
                    (b2 - a2) / offset,
                    linewidth=1.5,
                    edgecolor="green",
                    facecolor="green",
                    alpha=0.2,
                )
                ax.add_patch(rect_sym)"""

    if show_colorbar:
        if reserve_colorbar_space:
            fig.subplots_adjust(left=0.12, right=0.76, bottom=0.12, top=0.88)
        fig.colorbar(im, ax=ax, label=colorbar_label, pad=colorbar_pad)

    if legend_outside:
        # Keep overlay labels out of the data axes and reserve enough room for
        # all three handles without changing the matrix's square aspect.
        fig.tight_layout(rect=(0, 0, 0.76, 1))
    elif not reserve_colorbar_space:
        plt.tight_layout()
    if save_path:
        save_kwargs = {} if aspect == "equal" else {"bbox_inches": "tight"}
        fig.savefig(save_path, dpi=dpi, **save_kwargs)
        plt.close(fig)
    else:
        plt.show()


def read_bed_files(files: Union[str, List[str]]) -> List[pl.DataFrame]:
    if isinstance(files, str):
        files = [files]

    dataframes = []

    for file in files:
        with open(file, "r") as f:
            lines = [
                line.strip()
                for line in f
                if not line.startswith("track") and not line.startswith("#")
            ]

        if not lines:
            dataframes.append(pl.DataFrame())
            continue

        n_cols = len(lines[0].split("\t"))
        col_names = BED_COLUMNS[:n_cols] + [
            f"extra_{i}" for i in range(n_cols - len(BED_COLUMNS))
        ]

        rows = [line.split("\t")[:n_cols] for line in lines]
        df = pl.DataFrame(
            {col_names[i]: [row[i] for row in rows] for i in range(n_cols)}
        )

        # Cast specific columns to integer if present
        for col in ["start", "end", "thickStart", "thickEnd"]:
            if col in df.columns:
                df = df.with_columns(pl.col(col).cast(pl.Int64))

        dataframes.append(df)

    return dataframes


def top_n_frequent_distances(distances, n=5):
    # Count the occurrences of each distance
    distance_counts = Counter(distances)
    # Get the top n most common distances
    top_n = distance_counts.most_common(n)
    return top_n


def validate_json(path: str, required_keys: list = None) -> bool:
    if not os.path.isfile(path):
        print(f"[ERROR] File does not exist: {path}")
        return False

    try:
        with open(path, "r") as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        print(f"[ERROR] Invalid JSON syntax: {e}")
        return False

    if required_keys:
        missing = [k for k in required_keys if k not in data]
        if missing:
            print(f"[ERROR] Missing required keys in JSON: {missing}")
            return False

    return True


def write_summary_file(tuple_of_lists, out_csv_path: str) -> None:
    new_starts, new_ends, new_names, monomer, periodicity, hor = tuple_of_lists

    # Basic sanity check
    n = len(new_starts)
    if not (
        len(new_ends)
        == len(new_names)
        == len(monomer)
        == len(periodicity)
        == len(hor)
        == n
    ):
        raise ValueError("All lists in tuple_of_lists must have the same length.")

    # Group intervals by name, but treat None as unique per entry
    groups = {}  # key -> dict
    none_counter = 0

    for s, e, name, m, p, h in zip(
        new_starts, new_ends, new_names, monomer, periodicity, hor
    ):
        if name is None or name == "Unknown":
            none_counter += 1
            key = f"None_{none_counter}"  # unique row per None
            out_name = "Unclassified Repeat"
        else:
            key = name
            out_name = name

        if key not in groups:
            groups[key] = {
                "name": out_name,
                "monomer": m,
                "periodicity": p,
                "hor": h,
                "intervals": [],
            }

        groups[key]["intervals"].append((s, e))

    # Write CSV
    with open(out_csv_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["name", "monomer", "periodicity", "hor", "intervals"]
        )
        writer.writeheader()

        for _, row in groups.items():
            intervals_str = ";".join(f"{s}-{e}" for s, e in row["intervals"])
            writer.writerow(
                {
                    "name": row["name"],
                    "monomer": row["monomer"],
                    "periodicity": row["periodicity"],
                    "hor": row["hor"],
                    "intervals": intervals_str,
                }
            )
