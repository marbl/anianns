import argparse
from tokenize import group
from anianns.const import ASCII_ART, DESCRIPTION, VERSION
from itertools import islice
import polars as pl
import sys
import json
import os
import pysam
import math
import time
import numpy as np

from anianns.ani_matrix import (
    intersection_matrix,
    intersection_matrix_inverted,
    intersection_matrix_thresholded,
)

from anianns.build_kmer_db import (
    save_kmer_sets_shared_k,
)

from anianns.general_utils import (
    add_prefix_to_tuples,
    convert_dataframe_format,
    check_bed_vs_indexed_fasta,
    define_bounds,
    extract_region,
    extract_regions_by_name,
    extract_histograms_by_name,
    get_input_headers,
    plot_matrix,
    read_bed_files,
    validate_json,
    write_summary_file,
)

from anianns.kmer_utils import (
    build_kmer_sets,
    print_progress_bar,
)

from anianns.kmer_pipeline import SequenceBandPlan, iter_hashed_fasta_bands

from anianns.parse_matrix import (
    append_coordinates,
    get_diagonal_span,
    merge_shared_boundaries,
    sobel_with_diagonal_probes,
    sobel_spans,
    split_diagonal_attached,
)

from anianns.refine_boundaries import report_borders

from anianns.union_find import (
    sobel,
    sobel_with_diagonal_probes2,
    find_offdiag_rectangles,
)

os.environ["KMP_WARNINGS"] = "FALSE"


def mask_type(value):
    """Allow either integer or string values for --mask."""
    try:
        return int(value)
    except ValueError:
        return value


def validate_ntrprism_range(start, end, seq_len):
    """
    Validate the --range argument for the ntrprism subcommand.

    Returns an error message string on failure, or None on success.
    """
    if start >= end:
        return f"[ERROR] --range start ({start}) must be less than end ({end})."
    if start < 0 or end > seq_len:
        return f"[ERROR] --range ({start}, {end}) is out of bounds for sequence of length {seq_len}."
    return None


def format_matrix_runtime(index, total, runtime, previous_runtime=None):
    """Format a verbose per-matrix runtime and prior-matrix comparison."""
    prefix = f"Matrix {index}/{total} completed in {runtime:.3f} s"
    if previous_runtime is None:
        return f"{prefix} (comparison baseline)."

    difference = runtime - previous_runtime
    if abs(difference) < 0.0005:
        return f"{prefix} — same runtime as matrix {index - 1}."

    direction = "slower" if difference > 0 else "faster"
    percentage = (
        abs(difference) / previous_runtime * 100 if previous_runtime > 0 else 0.0
    )
    sign = "+" if difference > 0 else "-"
    return (
        f"{prefix} — {abs(difference):.3f} s {direction} "
        f"({sign}{percentage:.1f}%) than matrix {index - 1}."
    )


def get_parser():
    """
    Argument parsing for stand-alone runs.

    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=DESCRIPTION,
    )
    subparsers = parser.add_subparsers(
        dest="command", help="Choose mode: annotate, build_db, ntrprism"
    )
    annotate_parser = subparsers.add_parser(
        "annotate",
        help="Takes input fasta(s) and outputs an annotated bedfile of satellite arrays.",
    )
    build_db_parser = subparsers.add_parser(
        "build_db",
        help="Takes input fasta(s), bedfile(s) of known satellite coordinates, and a config file, and outputs a kme db (directory).",
    )
    ntrprism_parser = subparsers.add_parser(
        "ntrprism",
        help="Takes input fasta(s) plus regions of interest, outputs ntrprism k-mer spectra.",
    )
    annotate_parser.add_argument(
        "-f",
        "--fasta",
        default=argparse.SUPPRESS,
        help="Path to input fasta file(s).",
        required=True,
        nargs="+",
    )
    annotate_parser.add_argument(
        "-s",
        "--seq_id",
        nargs="+",
        default=None,
        help="Sequence ID to extract (multiple if using multifasta file). Will ignore if not found.",
    )
    annotate_parser.add_argument(
        "-d",
        "--directory",
        default=None,
        help="Name of output directory. Default: current working directory.",
    )
    annotate_parser.add_argument(
        "-o",
        "--output-format",
        choices=["bed", "gtf", "gff", "csv", "tsv", "json"],
        default="bed",
        help=(
            "Specify the output file format. "
            "Accepted values: bed, gtf, gff, csv, tsv, or json. "
            "Default: bed."
        ),
    )
    annotate_parser.add_argument(
        "-m",
        "--mask",
        nargs="*",
        type=mask_type,
        default=None,
        help=(
            "Name(s) or ID(s) of satellite arrays to mask. "
            "Use without parameters to apply default masking. "
            "Use 'ALL' for everything. Default: not applied."
        ),
    )
    annotate_parser.add_argument(
        "--soft",
        default=False,
        action="store_true",
        help="Apply soft masking (requires --mask).",
    )
    annotate_parser.add_argument(
        "-c",
        "--classify",
        default=None,
        help="Directory of satellite kmer db files, which Ani Ann's will use to classify into known satellite classes.",
    )
    annotate_parser.add_argument(
        "-t",
        "--threshold",
        default=None,
        help="Directory of satellite kmer db files, which Ani Ann's will use to classify into known satellite classes.",
    )
    annotate_parser.add_argument(
        "-k", "--kmer", type=int, default=21, help="k-mer length. Default: 21"
    )
    annotate_parser.add_argument(
        "--sketch",
        type=int,
        choices=(2, 4),
        default=4,
        help=(
            "Keep hashes divisible by this modulo when building window sketches. "
            "Use 2 for a denser, slower sketch or 4 for the default."
        ),
    )
    annotate_parser.add_argument(
        "-i", "--identity", type=int, default=86, help="Identity threshold. Default: 86"
    )
    annotate_parser.add_argument(
        "-w",
        "--window",
        type=int,
        default=2000,
        help="Dotplot window size, or the number of bp contained within each pixel in a plot. This is proportional to the sensitivity of satellite detection (ie. lower is more accurate, at the expense of runtime). Default: 2000.",
    )
    annotate_parser.add_argument(
        "--band",
        type=float,
        default=2.0,
        help="Max height in Mbp of band. Default: 2.0",
    )
    annotate_parser.add_argument(
        "--cache-dir",
        default=None,
        help=(
            "Directory for reusable canonical k-mer hash caches. "
            "Default: <output directory>/.anianns_cache."
        ),
    )
    annotate_parser.add_argument(
        "--identifier",
        help="Name of identifier. Used when no matches to a k-mer db are found, or if `--classify` is not provided. bed file to output to. Default: None",
    )
    annotate_parser.add_argument(
        "-p",
        "--plot",
        action="store_true",
        help="Output self-identity dotplot(s) for each sequence.",
    )
    annotate_group = annotate_parser.add_mutually_exclusive_group()
    annotate_group.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose logging output."
    )
    annotate_group.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress all logging output and text.",
    )

    build_db_parser.add_argument(
        "-f",
        "--fasta",
        default=argparse.SUPPRESS,
        required=True,
        help="Path to input fasta file(s).",
        nargs="+",
    )
    build_db_parser.add_argument(
        "-b",
        "--bed",
        default=argparse.SUPPRESS,
        required=True,
        help="Path to input bed file(s).",
        nargs="+",
    )
    build_db_parser.add_argument(
        "-c",
        "--config",
        help="Path to input config file(s).",
    )
    build_db_parser.add_argument(
        "-d",
        "--directory",
        default=None,
        help="Name of kmer db output directory. Default: In current working directory, will create a directory with the k-mer length.",
    )
    build_db_parser.add_argument(
        "-k", "--kmer", type=int, default=21, help="k-mer length. Default: 21"
    )
    build_db_group = build_db_parser.add_mutually_exclusive_group()
    build_db_group.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose logging output."
    )
    build_db_group.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress all logging output and text.",
    )

    ntrprism_parser.add_argument(
        "-f",
        "--fasta",
        default=argparse.SUPPRESS,
        help="Path to input fasta file(s).",
        required=True,
        nargs="+",
    )
    ntrprism_parser.add_argument(
        "-s",
        "--seq_id",
        nargs="+",
        default=None,
        help="Sequence ID to extract (multiple if using multifasta file). Will ignore if not found.",
    )
    ntrprism_parser.add_argument(
        "-r",
        "--range",
        default=None,
        type=int,
        nargs=2,
        metavar=("START", "END"),
        help="Genomic range as two integers (start end). If omitted, the full sequence is used.",
    )

    return parser


def main():
    print(ASCII_ART)
    print(f" {VERSION}")
    print("─" * 65)

    args = get_parser().parse_args()

    # -------- BUILD DB LOGIC --------#
    if args.command == "build_db":
        if not args.directory:
            foldername = f"k_{args.kmer}"
            args.directory = os.path.join(os.getcwd(), foldername)
        if not args.quiet:
            print(f"Building a k-mer database....\n")
        if args.config:
            with open(args.config, "r") as f:
                if not validate_json(args.config):
                    sys.exit(1)
                else:
                    try:
                        satellite_metadata = json.load(f)
                        print(f"Successfully loaded sat metadata from {args.config}!\n")
                    except:
                        print(
                            f"Error loading satellite metadata from  {args.config}. Exiting.\n"
                        )
                        sys.exit(1)
        else:
            satellite_metadata = None
            if not args.quiet:
                print(
                    f"No metadata provided. Creating a compressed kmer db for each unique satellite.\n"
                )

        if not args.quiet:
            print(f"Creating kmer db using {args.fasta}\n")
        bedfiles = read_bed_files(args.bed)

        check_bed_vs_indexed_fasta(bedfiles, args.fasta)

        # Normalize all BED DataFrames to have consistent schema before concatenating
        normalized_bedfiles = []
        required_columns = ["chrom", "start", "end", "name"]

        for df in bedfiles:
            # Ensure all required columns exist
            for col in required_columns:
                if col not in df.columns:
                    if col == "name":
                        # Use chrom as default name if missing
                        df = df.with_columns(pl.col("chrom").alias("name"))
                    else:
                        raise ValueError(
                            f"Required column '{col}' missing from BED file"
                        )

            # Select only the columns we need for consistency
            normalized_df = df.select(required_columns)
            normalized_bedfiles.append(normalized_df)

        # Concatenate normalized DataFrames
        bedfile_dfs = pl.concat(normalized_bedfiles)
        satellite_db = extract_regions_by_name(
            bedfile_dfs, args.fasta, args.kmer, args.verbose
        )
        # histogram_db = extract_histograms_by_name(bedfile_dfs,args.fasta,args.kmer,True)

        if satellite_metadata:
            for key in satellite_metadata:
                subtypes = [s.lower() for s in satellite_metadata[key]]
                print(subtypes)

                # Filter subtypes that exist in the database
                sub_kmer_db = {
                    name: satellite_db[name.lower()]
                    for name in subtypes
                    if name.lower() in satellite_db
                }

                if not sub_kmer_db:
                    continue  # Skip if none of the subtypes are found

                outputprefix = f"{key}.db"
                outputname = os.path.join(args.directory, outputprefix)
                for subtype_key, value in sub_kmer_db.items():
                    print(
                        f"Saving {subtype_key} into {outputname} containing {len(value)} k-mers"
                    )
                os.makedirs(args.directory, exist_ok=True)
                save_kmer_sets_shared_k(
                    sub_kmer_db, k=args.kmer, output_path=outputname
                )

        else:
            for key, value in satellite_db.items():
                subtypes = [key]
                sub_kmer_db = {name: satellite_db[name] for name in subtypes}
                outputprefix = f"{key}.db"
                outputname = os.path.join(args.directory, outputprefix)
                print(f"Saving {outputname} containing {len(value)} k-mers")
                os.makedirs(args.directory, exist_ok=True)
                save_kmer_sets_shared_k(
                    sub_kmer_db, k=args.kmer, output_path=outputname
                )

    # -------- ANNOTATE LOGIC --------#
    elif args.command == "annotate":
        missing_fastas = [fasta for fasta in args.fasta if not os.path.isfile(fasta)]
        if missing_fastas:
            for fasta in missing_fastas:
                print(f"[ERROR] FASTA file does not exist: {fasta}", file=sys.stderr)
            sys.exit(1)

        # Prep args
        band_height = int(args.band * 1_000_000)
        interval = (args.window + 1) // 2
        directory = args.directory or os.getcwd()
        hash_cache_dir = args.cache_dir or os.path.join(directory, ".anianns_cache")

        if not args.quiet:
            label_width = 20

            print(f"{'Output directory:':<{label_width}} {directory}")
            print(f"{'K-mer length:':<{label_width}} {args.kmer}")
            print(f"{'Sketch modulo:':<{label_width}} {args.sketch}")
            print(f"{'Band height:':<{label_width}} {band_height} bp")
            print(f"{'Window size:':<{label_width}} {args.window} bp")
            print(f"{'ANI threshold:':<{label_width}} {args.identity} %")
            print(f"{'Hash cache:':<{label_width}} {hash_cache_dir}")
            print(
                f"{'K-mer dir:':<{label_width}} {args.classify if args.classify else 'None provided'}"
            )
            print("─" * 65)

        # Build list of (fasta, [seq_ids]) pairs
        headers = get_input_headers(args.fasta)
        if args.seq_id:
            pairs = []
            for sid in args.seq_id:
                matches = [f for f, ids in headers if sid in ids]
                if matches:
                    pairs.append((matches[0], [sid]))
                else:
                    # Check if there's sequence bounds for this seq_id
                    seq_bounds = define_bounds(sid)
                    if seq_bounds:
                        print(seq_bounds)
                        matches = [f for f, ids in headers if seq_bounds[0] in ids]
                        if matches:
                            pairs.append((matches[0], [sid]))
                        else:
                            if not args.quiet:
                                print(f"Unable to locate {seq_bounds}. Skipping…")
                    else:
                        if not args.quiet:
                            print(f"Unable to locate {sid}. Skipping…")
        else:
            pairs = headers

        # 3) Open all FASTAs once
        fasta_handles = {f: pysam.FastaFile(f) for f, _ in pairs}
        try:
            # cache everything into locals
            k_param = args.kmer
            win = args.window
            verbosity = args.verbose
            build_sets = build_kmer_sets
            imat = intersection_matrix if args.plot else intersection_matrix_thresholded
            imat_inv = intersection_matrix_inverted
            get_span = get_diagonal_span
            merge_intv = merge_shared_boundaries

            # 4) Main loops
            for fasta, seq_ids in pairs:
                fh = fasta_handles[fasta]
                for seq_id in seq_ids:
                    # Define data structure for satellite coordinates. Variable names for sequence name, length, and if samtools was used for coordinates
                    satellite_coordinate_list = []
                    seq_len = fh.get_reference_length(seq_id)
                    seq_bounds = define_bounds(seq_id)

                    if seq_bounds and not args.quiet:
                        print(f"Found bounds for {seq_id}: {seq_bounds}\n")

                    print(f"Creating an ANI matrix for {seq_id}:\n")

                    if seq_len < k_param:
                        if not args.quiet:
                            print(
                                f"Sequence {seq_id} is shorter than k={k_param}; skipping.\n"
                            )
                        continue

                    # If seq_len is close to the band height, then keep everything in the same band.
                    sequence_band_height = band_height
                    if seq_len < band_height + (band_height // 2):
                        sequence_band_height = seq_len
                        if not args.quiet:
                            print(f"Adjusting band to {sequence_band_height} bp.\n")

                    band_plan = SequenceBandPlan(
                        sequence_length=seq_len,
                        kmer=k_param,
                        band_height=sequence_band_height,
                        windows=(win,),
                    )
                    n_windows = band_plan.band_count
                    max_len = band_plan.window_plans[win].max_len
                    band_iterator = iter(
                        iter_hashed_fasta_bands(
                            fasta,
                            seq_id,
                            band_plan,
                            cache_dir=hash_cache_dir,
                        )
                    )
                    try:
                        first_band = next(band_iterator)
                    except StopIteration:
                        continue

                    if not args.quiet:
                        print_progress_bar(
                            0,
                            n_windows,
                            prefix="Progress:",
                            suffix="Complete",
                            length=40,
                        )

                    # Create initial window
                    kmers_list = band_plan.hashes_for_window(first_band, win)
                    prev_ov, prev_nov = build_sets(
                        kmers_list, max_len, win, interval, sketch=args.sketch
                    )
                    matrix_started = time.perf_counter()
                    initial_matrix = (
                        imat(prev_ov, prev_nov, k_param)
                        if args.plot
                        else imat(prev_ov, prev_nov, k_param, args.identity)
                    )
                    previous_matrix_runtime = time.perf_counter() - matrix_started
                    if verbosity:
                        print(
                            format_matrix_runtime(
                                1, n_windows, previous_matrix_runtime
                            )
                        )

                    if n_windows > 1:
                        if not args.quiet:
                            print_progress_bar(
                                1,
                                n_windows,
                                prefix="Progress:",
                                suffix="Complete",
                                length=40,
                            )
                        # Get spans for the initial window
                        initial_detection_matrix = (
                            initial_matrix >= args.identity
                            if args.plot
                            else initial_matrix
                        )
                        spans = get_span(initial_detection_matrix, win, zero_tol=2)

                        # TODO: Remove low count spans
                        """for element in spans:
                            print(element, element[1], element[1]*win, element[0][1]-element[0][0])"""

                        # Replace 0 here with start prefix
                        for coordinates in merge_intv(
                            intervals=spans, prefix=0, window=win, verbose=False
                        ):
                            satellite_coordinate_list.append(coordinates)
                        # Iterate through independently fetched/hash-overlapped bands.
                        # The producer queues one future band while this process
                        # performs the current band's matrix work.
                        for hashed_band in band_iterator:
                            w = hashed_band.index + 1
                            kmers_list = band_plan.hashes_for_window(hashed_band, win)

                            ov, nov = build_sets(
                                kmers_list, max_len, win, interval, sketch=args.sketch
                            )

                            matrix_started = time.perf_counter()
                            updated_matrix = (
                                imat(ov, nov, k_param)
                                if args.plot
                                else imat(ov, nov, k_param, args.identity)
                            )
                            matrix_runtime = time.perf_counter() - matrix_started
                            if verbosity:
                                print(
                                    format_matrix_runtime(
                                        w,
                                        n_windows,
                                        matrix_runtime,
                                        previous_matrix_runtime,
                                    )
                                )
                            previous_matrix_runtime = matrix_runtime

                            if args.plot:
                                inv = imat_inv(
                                    initial_matrix,
                                    updated_matrix,
                                    prev_ov,
                                    prev_nov,
                                    ov,
                                    nov,
                                    k_param,
                                )
                                inv[inv < args.identity] = 0
                                plot_matrix(inv)
                                os.makedirs(directory, exist_ok=True)
                                matrix_filename = os.path.join(
                                    directory, f"{seq_id}_{w}_matrix.npy"
                                )
                                np.save(matrix_filename, inv)
                                if not args.quiet:
                                    print(f"Saved matrix to {matrix_filename}")
                                M_diag, M_distal = split_diagonal_attached(inv)
                                diag_filename = os.path.join(
                                    directory, f"{seq_id}_diag_{w}_matrix.npy"
                                )
                                np.save(diag_filename, M_diag)
                                distal_filename = os.path.join(
                                    directory, f"{seq_id}_distal_{w}_matrix.npy"
                                )
                                np.save(distal_filename, M_distal)

                            updated_detection_matrix = (
                                updated_matrix >= args.identity
                                if args.plot
                                else updated_matrix
                            )
                            new_spans = get_span(
                                updated_detection_matrix, win, zero_tol=2
                            )
                            prefix_amount = hashed_band.start

                            # print(new_spans)
                            """if verbosity:
                                print(f"Current prefix: {prefix_amount}\n")"""
                            for coordinates in merge_intv(
                                intervals=new_spans,
                                prefix=prefix_amount,
                                window=win,
                                verbose=False,
                            ):
                                satellite_coordinate_list.append(coordinates)
                            """probes = sobel_with_diagonal_probes(
                                M = inv,
                                thresh = 0.7,
                                min_thick = 1,
                            )"""
                            # print(sobel_spans(probes,win))
                            # Roll matrices forward
                            initial_matrix, prev_ov, prev_nov = updated_matrix, ov, nov

                            # Update progress bar
                            if not args.quiet:
                                if w == n_windows:
                                    print_progress_bar(
                                        n_windows,
                                        n_windows,
                                        prefix="Progress:",
                                        suffix="Completed!\n",
                                        length=40,
                                    )
                                else:
                                    print_progress_bar(
                                        w,
                                        n_windows,
                                        prefix="Progress:",
                                        suffix="Complete",
                                        length=40,
                                    )

                    else:
                        # No progress bar in this case
                        initial_detection_matrix = (
                            initial_matrix >= args.identity
                            if args.plot
                            else initial_matrix
                        )
                        spans = get_span(initial_detection_matrix, win, zero_tol=2)

                        if args.plot:
                            initial_matrix[initial_matrix < args.identity] = 0
                            plot_matrix(initial_matrix)
                            M_diag, M_distal = split_diagonal_attached(initial_matrix)

                        for coordinates in merge_intv(
                            intervals=spans, prefix=0, window=win, verbose=False
                        ):
                            # Append only if the counts pass a threshold
                            # Counting threshold function here
                            satellite_coordinate_list.append(coordinates)

                    filtered = []
                    for x, y, count in satellite_coordinate_list:
                        size = y - x
                        count_size = count * win
                        if (count_size <= size * 0.6) or (
                            count < 10 and count_size <= size * 0.75
                        ):
                            continue
                        filtered.append((x, y, count))
                    satellite_coordinate_list = filtered
                    # print(satellite_coordinate_list)
                    # sys.exit(0)
                    if seq_bounds:
                        # seq_bounds[0] is the name seq_bounds[1] is the start offset, 2 is the end offset
                        merged_coordinates = [
                            (
                                seq_bounds[0],  # chrom
                                start + int(seq_bounds[1]),
                                end + int(seq_bounds[1]) - 1,
                                seq_bounds[0],  # name
                                0,  # score
                                ".",  # strand
                                start + int(seq_bounds[1]),  # thickStart
                                # end + seq_bounds[1],      # thickEnd
                                end + int(seq_bounds[1]) - 1,  # thickEnd
                                "255,0,0",  # itemRgb
                            )
                            for start, end, _ in filtered
                        ]
                    else:
                        merged_coordinates = [
                            (
                                seq_id,  # chrom
                                start,  # start
                                end,  # end
                                seq_id,  # name
                                0,  # score
                                ".",  # strand
                                start,  # thickStart
                                end,  # thickEnd
                                "255,0,0",  # itemRgb
                            )
                            for start, end, _ in filtered
                        ]
                    """if verbosity:
                        print(merged_coordinates)"""

                    df1 = pl.DataFrame(
                        merged_coordinates,
                        schema=[
                            "#chrom",
                            "start",
                            "end",
                            "name",
                            "score",
                            "strand",
                            "thickStart",
                            "thickEnd",
                            "itemRgb",
                        ],
                        orient="row",
                    )

                    # TODO: Fix formatting
                    '''if args.output_format == "bed":
                        suffix = "bed"
                    if args.output_format != "bed":
                        if args.output_format == "gtf":
                            df_converted = convert_dataframe_format(df1, "gtf")
                            suffix = "gtf"
                        elif args.output_format == "gff":
                            df_converted = convert_dataframe_format(df1, "gff")
                            suffix = "gff"
                        elif args.output_format == "csv":
                            df_converted = convert_dataframe_format(df1, "csv")
                            suffix = "csv"
                        elif args.output_format == "tsv":
                            df_converted = convert_dataframe_format(df1, "tsv")
                            suffix = "tsv"
                        elif args.output_format == "json":
                            df_converted = convert_dataframe_format(df1, "json")
                            suffix = "json"
                        else:
                            sys.exit(f"[ERROR] Unknown output format: {args.output_format}. Defaulting to bed.\n")
                            suffix = "bed"
                    else:
                        suffix = "bed"'''
                    suffix = "bed"
                    df_converted = df1
                    """annotation_file_name = f"{seq_id}_unrefined.{suffix}"
                    annotation_file_path = os.path.join(directory, annotation_file_name)
                    os.makedirs(directory, exist_ok=True)
                    df_converted.write_csv(annotation_file_path, separator="\t")"""

                    # If we are using a subseqeunce of a larger fasta, we need to adjust the coordinates back to the original reference frame before outputting
                    if seq_bounds:
                        # fa, seq_id, seq_len, band, offset, window, k, df: pl.DataFrame, classify, verbose, quiet) -> None:
                        tuple_of_lists = report_borders(
                            fa=fh,
                            seq_id=seq_id,
                            seq_len=seq_len,
                            band=args.band,
                            offset=seq_bounds[1],
                            window=win,
                            k=k_param,
                            df=df1,
                            classify=args.classify,
                            verbose=verbosity,
                            quiet=args.quiet,
                        )
                        (
                            new_starts,
                            new_ends,
                            new_names,
                            monomer,
                            periodicity,
                            hor,
                        ) = tuple_of_lists
                        # print(tuple_of_lists)
                    else:
                        tuple_of_lists = report_borders(
                            fa=fh,
                            seq_id=seq_id,
                            seq_len=seq_len,
                            band=args.band,
                            offset=1,
                            window=win,
                            k=k_param,
                            df=df1,
                            classify=args.classify,
                            verbose=verbosity,
                            quiet=args.quiet,
                        )
                        (
                            new_starts,
                            new_ends,
                            new_names,
                            monomer,
                            periodicity,
                            hor,
                        ) = tuple_of_lists
                        # print(tuple_of_lists)

                    # Replace the columns in df1
                    if seq_bounds:
                        offset = int(seq_bounds[1])
                        df2 = pl.DataFrame(
                            {
                                "#chrom": [seq_id] * len(new_starts),
                                "start": [s + offset for s in new_starts],
                                "end": [e + offset for e in new_ends],
                                "name": new_names,
                                "score": [e for e in monomer],
                                "strand": ["."] * len(new_starts),
                                "thickStart": [s + offset for s in new_starts],
                                "thickEnd": [e + offset for e in new_ends],
                                "itemRgb": ["0,0,0"] * len(new_starts),
                            }
                        )
                    else:
                        df2 = pl.DataFrame(
                            {
                                "#chrom": [seq_id] * len(new_starts),
                                "start": new_starts,
                                "end": new_ends,
                                "name": new_names,
                                "score": [e for e in monomer],
                                "strand": ["."] * len(new_starts),
                                "thickStart": new_starts,
                                "thickEnd": new_ends,
                                "itemRgb": ["0,0,0"] * len(new_starts),
                            }
                        )

                    bedfilename = f"{seq_id}.bed"
                    bedfilepath = os.path.join(directory, bedfilename)
                    os.makedirs(directory, exist_ok=True)
                    df2.write_csv(bedfilepath, separator="\t")

                    csvfilename = f"{seq_id}.csv"
                    csvfilepath = os.path.join(directory, csvfilename)
                    write_summary_file(tuple_of_lists, csvfilepath)

                    print(
                        f"Successfully finished annotating {seq_id} to {bedfilepath}\n"
                    )

        except Exception as e:
            print(e)

    # -------- NTRPRISM LOGIC --------#
    elif args.command == "ntrprism":
        if args.range is not None:
            start, end = args.range
            for fasta_path in args.fasta:
                fh = pysam.FastaFile(fasta_path)
                seq_ids = args.seq_id if args.seq_id else fh.references
                for seq_id in seq_ids:
                    seq_len = fh.get_reference_length(seq_id)
                    error = validate_ntrprism_range(start, end, seq_len)
                    if error:
                        print(error)
                        sys.exit(1)
                fh.close()

            headers = get_input_headers(args.fasta)
        if args.seq_id:
            pairs = []
            for sid in args.seq_id:
                matches = [f for f, ids in headers if sid in ids]
                if matches:
                    pairs.append((matches[0], [sid]))
                else:
                    # Check if there's sequence bounds for this seq_id
                    seq_bounds = define_bounds(sid)
                    if seq_bounds:
                        # print(seq_bounds)
                        matches = [f for f, ids in headers if seq_bounds[0] in ids]
                        if matches:
                            pairs.append((matches[0], [sid]))
                        else:
                            if not args.quiet:
                                print(f"Unable to locate {seq_bounds}. Skipping…")
                    else:
                        if not args.quiet:
                            print(f"Unable to locate {sid}. Skipping…")
        else:
            pairs = headers

        # 3) Open all FASTAs once
        fasta_handles = {f: pysam.FastaFile(f) for f, _ in pairs}
