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

from anianns.ani_matrix import (
    intersection_matrix,
    intersection_matrix_inverted
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
)

from anianns.kmer_utils import (
    build_kmer_sets,
    generate_kmers_from_fasta,
    print_progress_bar
)

from anianns.parse_matrix import (
    append_coordinates,
    get_diagonal_span,
    merge_shared_boundaries,
    split_diagonal_attached
)

from anianns.refine_boundaries import report_borders

os.environ["KMP_WARNINGS"] = "FALSE"

def mask_type(value):
    """Allow either integer or string values for --mask."""
    try:
        return int(value)
    except ValueError:
        return value


def get_parser():
    """
    Argument parsing for stand-alone runs.

    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=DESCRIPTION,
    )
    subparsers = parser.add_subparsers(
        dest="command", help="Choose mode: annotate or mask"
    )
    annotate_parser = subparsers.add_parser(
        "annotate",
        help="Takes input fasta(s) and outputs an annotated bedfile of satellite arrays.",
    )
    build_db_parser = subparsers.add_parser(
        "build_db",
        help="Takes input fasta(s), bedfile(s) of known satellite coordinates, and a config file, and outputs a kme db (directory).",
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
        "-o", "--output-format",
        choices=["bed", "gtf", "gff", "csv", "tsv", "json"],
        default="bed",
        help=(
            "Specify the output file format. "
            "Accepted values: bed, gtf, gff, csv, tsv, or json. "
            "Default: bed."
        ),
    )
    annotate_parser.add_argument(
        "-m", "--mask",
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
        "-v", "--verbose", action="store_true",
        help="Enable verbose logging output."
    )
    annotate_group.add_argument(
        "-q", "--quiet", action="store_true",
        help="Suppress all logging output and text."
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
        "-v", "--verbose", action="store_true",
        help="Enable verbose logging output."
    )
    build_db_group.add_argument(
        "-q", "--quiet", action="store_true",
        help="Suppress all logging output and text."
    )

    return parser


def main():
    print(ASCII_ART)
    print(f" {VERSION}")
    print("─" * 65)

    args = get_parser().parse_args()

    #-------- BUILD DB LOGIC --------#
    if args.command == "build_db":
        if not args.directory:
            foldername = f"k_{args.kmer}"
            args.directory = os.path.join(os.getcwd(),foldername)
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
                        raise ValueError(f"Required column '{col}' missing from BED file")
            
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

    #-------- ANNOTATE LOGIC --------#
    elif args.command == "annotate":

        # Prep args
        band_height = int(args.band * 1_000_000)
        interval    = (args.window + 1) // 2
        directory     = args.directory or os.getcwd()

        if not args.quiet:
            print(f"Output directory: {directory}")
            print(f"K-mer length:     {args.kmer}")
            print(f"Band height:      {band_height} bp")
            print(f"Window size:      {args.window} bp")
            print(f"ANI threshold:    {args.identity} %")
            if args.classify:
                print(f"K-mer dir:   {args.classify}")
            else:
                print(f"K-mer dir:   None provided")
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
            k_param    = args.kmer
            win        = args.window
            build_sets = build_kmer_sets
            imat       = intersection_matrix
            imat_inv   = intersection_matrix_inverted
            get_span   = get_diagonal_span
            merge_intv = merge_shared_boundaries

            window_size = band_height + interval
            max_len     = (band_height + win) // win

            # 4) Main loops
            for fasta, seq_ids in pairs:
                fh = fasta_handles[fasta]
                for seq_id in seq_ids:
                    # Define data structure for satellite coordinates. Variable names for sequence name, length, and if samtools was used for coordinates
                    satellite_coordinate_list = []
                    seq_str = fh.fetch(seq_id)
                    seq_len = len(seq_str)
                    seq_bounds = define_bounds(seq_id)

                    if seq_bounds and not args.quiet:
                        print(f"Found bounds for {seq_id}: {seq_bounds}\n")

                    print(f"Creating an ANI matrix for {seq_id}:\n")

                    # If seq_len is close to the band height, then keep everything in the same band.
                    if seq_len < band_height + (band_height // 2):
                        band_height = seq_len
                        n_windows = 1
                        if seq_len > band_height:
                            print(f"Adjusting band to {band_height} bp.\n")
                    else:
                        n_windows = math.ceil(seq_len / band_height)
                        remainder = seq_len % n_windows
                        if remainder > 0 and remainder < band_height // 10:
                            n_windows -= 1

                    # k‑mer iterator. This is the generator function.
                    kmer_it = generate_kmers_from_fasta(seq_str, k_param, True)

                    if not args.quiet:
                        print_progress_bar(0, n_windows, prefix="Progress:", suffix="Complete", length=40)

                    # Create initial window
                    kmers_list = list(islice(kmer_it, window_size))
                    prev_ov, prev_nov = build_sets(kmers_list, max_len, win, interval)
                    initial_matrix = imat(prev_ov, prev_nov, k_param)

                    # Split matrix into M_diag and M_distal
                    M_diag, M_distal = split_diagonal_attached(initial_matrix)

                    # Plot here
                    if args.verbose:
                        print(f"Partitioning into {n_windows} windows of {band_height} bp each.\n")
                    if n_windows > 1:
                        if not args.quiet:
                            print_progress_bar(1, n_windows, prefix="Progress:", suffix="Complete", length=40)
                        # Get spans for the initial window
                        spans = get_span(initial_matrix, win, zero_tol=2)

                        # TODO: Remove low count spans
                        '''for element in spans:
                            print(element, element[1], element[1]*win, element[0][1]-element[0][0])'''
                        
                        # Replace 0 here with start prefix
                        for coordinates in merge_intv(spans, 0, False):
                                satellite_coordinate_list.append(coordinates)
                        # Iterate through the remaining windows
                        for w in range(2, n_windows+1):
                            if w == n_windows:
                                chunk = list(kmer_it)
                            else:
                                chunk = list(islice(kmer_it, band_height))
                            if not chunk:
                                break

                            # Keep last interval from old + new chunk
                            kmers_list = kmers_list[-interval:] + chunk

                            ov, nov = build_sets(kmers_list, max_len, win, interval)
                            
                            updated_matrix = imat(ov, nov, k_param)

                            inv = imat_inv(
                                initial_matrix, updated_matrix,
                                prev_ov, prev_nov,
                                ov, nov,
                                k_param
                            )
                            inv[inv < args.identity] = 0

                            '''plot_matrix(inv)'''

                            new_spans = get_span(updated_matrix, win, zero_tol=2)
                            prefix_amount = band_height * (w-1)
                            
                            '''if args.verbose:
                                print(f"Current prefix: {prefix_amount}\n")'''

                            for coordinates in merge_intv(new_spans, prefix_amount, False):
                                satellite_coordinate_list.append(coordinates)

                            # Roll matrices forward
                            initial_matrix, prev_ov, prev_nov = updated_matrix, ov, nov

                            # Update progress bar
                            if not args.quiet:
                                if w == n_windows:
                                    print_progress_bar(n_windows, n_windows, prefix="Progress:", suffix="Completed!\n", length=40)
                                else:
                                    print_progress_bar(w, n_windows, prefix="Progress:", suffix="Complete", length=40)

                            if args.verbose:
                                print(new_spans)

                    else:
                        # No progress bar in this case
                        initial_matrix[initial_matrix < args.identity] = 0
                        spans = get_span(initial_matrix, win, zero_tol=2)

                        # Get off diagonals here after screening for identity
                        M_diag, M_offdiag = split_diagonal_attached(initial_matrix)
                        for coordinates in merge_intv(spans, 0, False):
                            satellite_coordinate_list.append(coordinates)

                    filtered_spans = []
                    seen_x = set()
                    seen_y = set()
                    for x, y in satellite_coordinate_list:
                        if x in seen_x or y in seen_y:
                            continue
                        seen_x.add(x)
                        seen_y.add(y)
                        filtered_spans.append((x, y))

                    '''if args.verbose:
                        print(f"{filtered_spans} filtered spans\n")'''

                    if seq_bounds:
                        # seq_bounds[0] is the name seq_bounds[1] is the start offset, 2 is the end offset
                        merged_coordinates = [
                            (
                                seq_bounds[0],   # chrom
                                start + int(seq_bounds[1]),
                                end + int(seq_bounds[1]) - 1,
                                seq_bounds[0],   # name
                                0,        # score
                                ".",      # strand
                                start + int(seq_bounds[1]),    # thickStart
                                #end + seq_bounds[1],      # thickEnd
                                end + int(seq_bounds[1]) - 1, # thickEnd
                                "255,0,0"       # itemRgb
                            )
                            for start, end in filtered_spans
                        ]
                    else:
                        merged_coordinates = [
                            (
                                seq_id,   # chrom
                                start,    # start
                                end,      # end
                                seq_id,   # name
                                0,        # score
                                ".",      # strand
                                start,    # thickStart
                                end,      # thickEnd
                                "255,0,0" # itemRgb
                            )
                            for start, end in filtered_spans
                        ]
                    '''if args.verbose:
                        print(merged_coordinates)'''

                    df1 = pl.DataFrame(
                        merged_coordinates,
                        schema=[
                            "#chrom", "start", "end",
                            "name", "score", "strand",
                            "thickStart", "thickEnd", "itemRgb"
                        ],
                        orient="row"
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
                    '''annotation_file_name = f"{seq_id}_unrefined.{suffix}"
                    annotation_file_path = os.path.join(directory, annotation_file_name)
                    os.makedirs(directory, exist_ok=True)
                    df_converted.write_csv(annotation_file_path, separator="\t")'''

                    if seq_bounds:
                        #fa, seq_id, seq_len, band, offset, window, k, df: pl.DataFrame, classify, verbose, quiet) -> None:
                        tuple_of_lists = report_borders(
                            fa=fasta,
                            seq_id=seq_id,
                            seq_len=seq_len, 
                            band=args.band,
                            offset=seq_bounds[1],
                            window=win,
                            k=k_param,
                            df=df1,
                            classify=args.classify,
                            verbose=args.verbose,
                            quiet=args.quiet
                        )
                        new_starts, new_ends, new_names = tuple_of_lists
                        if args.verbose:
                            print(tuple_of_lists)
                    else:
                        tuple_of_lists = report_borders(
                            fa=fasta,
                            seq_id=seq_id,
                            seq_len=seq_len, 
                            band=args.band,
                            offset=1,
                            window=win,
                            k=k_param,
                            df=df1,
                            classify=args.classify,
                            verbose=args.verbose,
                            quiet=args.quiet
                        )
                        new_starts, new_ends, new_names = tuple_of_lists
                        if args.verbose:
                            print(tuple_of_lists)

                    # Replace the columns in df1
                    df2 = pl.DataFrame({
                        "chrom": [seq_id] * len(new_starts),
                        "start": new_starts,
                        "end": new_ends,
                        "name": new_names,
                        "score": [0] * len(new_starts),
                        "strand": ["."] * len(new_starts),
                        "thickStart": new_starts,
                        "thickEnd": new_ends,
                        "itemRgb": ["0,0,0"] * len(new_starts)
                    })

                    bedfilename = f"{seq_id}.bed"
                    bedfilepath = os.path.join(directory,bedfilename)
                    os.makedirs(directory,exist_ok=True)
                    df2.write_csv(bedfilepath, separator="\t")

                    print(f"Successfully finished annotating {seq_id} to {bedfilepath}\n")

        except Exception as e:
            print(e)
