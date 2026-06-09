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
import csv

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

    # --- BED → GTF conversion ---
    if format == "gtf":
        # GTF columns: seqname, source, feature, start, end, score, strand, frame, attribute
        gtf_df = df.select(
            [
                pl.col("chrom").alias("seqname"),
                pl.lit("converted").alias("source"),
                pl.lit("exon").alias("feature"),
                pl.col("start"),
                pl.col("end"),
                pl.col("score"),
                pl.col("strand"),
                pl.lit(".").alias("frame"),
                (
                    pl.concat_str([pl.lit('gene_id "'), pl.col("name"), pl.lit('";')])
                ).alias("attribute"),
            ]
        )
        return gtf_df.write_csv(None, separator="\t", include_header=False)

    # --- BED → GFF conversion ---
    elif format == "gff":
        # GFF columns: seqid, source, type, start, end, score, strand, phase, attributes
        gff_df = df.select(
            [
                pl.col("chrom").alias("seqid"),
                pl.lit("converted").alias("source"),
                pl.lit("region").alias("type"),
                pl.col("start"),
                pl.col("end"),
                pl.col("score"),
                pl.col("strand"),
                pl.lit(".").alias("phase"),
                (pl.concat_str([pl.lit("ID="), pl.col("name")])).alias("attributes"),
            ]
        )
        return gff_df.write_csv(None, separator="\t", include_header=False)

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
        fasta_file (str): Path to the FASTA file (indexed with .fai).
        chr (str): Chromosome or contig name.
        region_start (int): 1-based start coordinate.
        region_end (int): 1-based end coordinate (inclusive).

    Returns:
        str or None: Extracted DNA sequence, or None if an error occurred.
    """
    try:
        if region_start < 1:
            region_start = 1
        fasta = pysam.FastaFile(fasta_file)
        sequence = fasta.fetch(chr, region_start, region_end)
        fasta.close()
        return sequence
    except Exception as e:
        print(
            f"Error fetching region {chr}:{region_start}-{region_end} from {fasta_file}\n"
            f"Details: {e}\n"
        )
        return None


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
):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    if not isinstance(matrix, np.ndarray):
        raise TypeError("Input must be a NumPy array")
    if matrix.ndim != 2:
        raise ValueError("Input must be a 2D matrix")

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    im = ax.imshow(matrix, cmap=cmap, aspect="auto", vmin=86, interpolation="nearest")

    ax.set_title(title)
    ax.set_xlabel("Columns")
    ax.set_ylabel("Rows")

    # --- scale tick labels only ---
    xticks = ax.get_xticks()
    yticks = ax.get_yticks()

    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels([f"{x * offset:.2f}" for x in xticks])
    ax.set_yticklabels([f"{y * offset:.2f}" for y in yticks])

    # --- multiple highlight regions ---
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
            )
            ax.add_patch(rect)

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
        plt.colorbar(im, ax=ax, label="Value")

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")

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
