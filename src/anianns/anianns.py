#!/usr/bin/env python3
import argparse
import math
import os
import numpy as np
import sys
from moddotplot.parse_fasta import readKmersFromFile, getInputHeaders, isValidFasta
from moddotplot.estimate_identity import createSelfMatrix
from anianns.parse_input import extract_region, read_bedfile, readKmersFromFileNonHashed, readSequenceKmersFromFile, readSequenceKmersFromFileForward, generateKmersFromFastaForward, generateKmersFromFastaReverse
from anianns.parse_dotplot import (
    mask_fasta_with_bed,
    row_non_zero_stats,
    save_to_file,
    find_non_zero_length,
    mask_offdiagonals,
    mask_diagonals,
    merge_coordinates,
    create_bed_file,
)
from anianns.parse_bedfile import classify_seq, get_color, detect_transitions, determine_strand, find_files_with_suffix, get_unique_first_column_values
from collections import Counter
from datetime import datetime
from itertools import islice
from anianns.const import COLOR_MAPPING


def get_parser():
    """
    Argument parsing for stand-alone runs.

    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Ani Ann's: Ani augmented Annotation of satellite arrays",
    )
    parser.add_argument(
        "-f",
        "--fasta",
        default=argparse.SUPPRESS,
        help="Path to input fasta file(s).",
        nargs="+",
    )
    parser.add_argument(
        "--band", type=float, default=8.0, help="Max height in Mbp of band."
    )
    parser.add_argument("-k", "--kmer", default=21, help="k-mer length")
    parser.add_argument(
        "--overlap", type=float, default=0.1, help="Percent overlap. Must be < 0.5."
    )
    parser.add_argument("--identity", type=int, default=86, help="Identity threshold.")
    parser.add_argument(
        "-o",
        "--output-dir",
        default=None,
        help="Directory name for saving matrices and coordinate logs. Defaults to working directory.",
    )
    parser.add_argument(
        "-w", "--window", type=int, default=2000, help="Window size of ModDotPlot."
    )
    parser.add_argument(
        "-m", "--mask", action="store_true", help="Create a masked fasta file."
    )
    parser.add_argument(
        "-c", "--classify", help="Classify satellites into known satellite classes."
    )
    parser.add_argument(
        "-b",
        "--bed",
        default=None,
        help="Name of output bed file. Default is anianns_YYYYMMDD.bed",
    )
    return parser


def main():
    args = get_parser().parse_args()
    skip_fasta = False
    # Check for command line arg errors:
    if not args.fasta and not args.bed:
        print("Either -f/--fasta or -b/--bed must be provided.")
        sys.exit(1)

    if args.bed:
        skip_fasta = True

    ASCII_ART = """
     _          _      _                    
    / \   _ __ (_)    / \   _ __  _ __  ' ___ 
   / _ \ | '_ \| |   / _ \ | '_ \| '_ \  / __|
  / ___ \| | | | |  / ___ \| | | | | | | \__ |
 /_/   \_\_| |_|_| /_/   \_\_| |_|_| |_| |___/
    """
    print(ASCII_ART)
    if skip_fasta:
        if args.classify:
            seq_bed_list = []
            fasta_list = args.fasta.copy()
            print("Classifying")
            for i in args.fasta:
                try:
                    isInBed = False
                    isValidFasta(i)
                    headers = getInputHeaders(i)
                    # Check that fasta file and bed file matchup so that we're not wasting time.
                    for head in headers:
                        print(head)
                    names = get_unique_first_column_values(args.bed)
                    print(names)
                    for k in names:
                        if k in headers:
                            print(f"{k} found in fasta file {i}")
                            isInBed = True
                    if isInBed:
                        seq_bed_list.append((headers, i, args.bed))
                except Exception as e:
                    print(
                        f"\nUnable to open {i}. Please check it is correctly formatted or compressed...\n"
                    )
                fasta_list.remove(i)
            sat_db = find_files_with_suffix(args.classify,".bin")
            print("Finsihed sats")
            if len(seq_bed_list) < 1:
                print(f"Not found in fasta file. Double check that names match up.")
            for thing in seq_bed_list:
                print(thing[1],thing[2])
                classify_seq(thing[1],thing[2],args.kmer,sat_db)
            sys.exit(0)
        # Make sure the matching fasta file is included in order to get orientation
        print(args.fasta[0])
        print(args.bed)
        print(determine_strand(args.bed, args.window, args.fasta[0], args.kmer))

    else:   
        # -----------FASTA INPUT SEQUENCE VALIDATION-----------
        seq_list = []
        fasta_list = args.fasta.copy()
        for i in args.fasta:
            try:
                isValidFasta(i)
                headers = getInputHeaders(i)

                if len(headers) > 1:
                    print(f"File {i} contains multiple fasta entries.\n")

                seq_list.append((headers, i))  # Add all headers to seq_list

            except Exception as e:
                print(
                    f"\nUnable to open {i}. Please check it is correctly formatted or compressed...\n"
                )
                fasta_list.remove(i)

        unique_sequences = [x[0] for x in seq_list]
        flattened_unique_sequences = [
            item for sublist in unique_sequences for item in sublist
        ]

        if len(flattened_unique_sequences) > 1:
            print(f"Annotating the following sequences: {flattened_unique_sequences}. \n")
        else:
            print(f"Annotating {flattened_unique_sequences}: \n")

        # -----------CREATE PATH IF NOT GIVEN------------
        if not args.output_dir:
            args.output_dir = os.getcwd()
        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)

        # -----------Run Ani Ann's on sequences-----------

        if not args.bed:
            timestamp = datetime.now().strftime("%Y%m%d")
            bed_name = f"anianns_{timestamp}.bed"
        else:
            bed_name = args.bed
        bed_coordinates = []
        print(f"Writing annotations to {os.path.join(args.output_dir,bed_name)}:\n")

        for fasta_file in seq_list:
            fasta_name = fasta_file[1]
            for sequence in fasta_file[0]:
                # Initialize bed file db. i[1] = fa file, i[0] = seq header
                k = readSequenceKmersFromFile(fasta_name, sequence, args.kmer, False)
                genome_length = len(k[0])
                divisor = math.ceil(genome_length / (args.band * 1000000))
                for j in range(0, divisor):
                    start = int(j * 1000000 * args.band)
                    end = int(min(start + (1000000 * args.band), genome_length))
                    print(
                        f"Annotating {sequence} from {start} to {end} ({j+1}/{divisor}): \n"
                    )
                    submatrix = createSelfMatrix(
                        sequence_length = len(k[0][start:end]),
                        sequence = k[0][start:end],
                        window_size = args.window,
                        sparsity = 4,
                        delta = 0.5,
                        k = args.kmer,
                        identity = args.identity,
                        ambiguous = False,
                        sketch_size = 500,
                    )
                    nonzerodata = find_non_zero_length(submatrix, args.window)

                    tuple_counts = Counter(nonzerodata[1])
                    merged_result = merge_coordinates(tuple_counts)
                    bed_coordinates = []
                    for (x, y), count in merged_result:
                        x_coord = int((j * 1000000 * args.band) + x)
                        y_coord = int((j * 1000000 * args.band) + y)
                        bed_coordinates.append((sequence, x_coord, y_coord, 'censat', 0, ".", x_coord, y_coord, COLOR_MAPPING["other"]))
                    create_bed_file(bed_coordinates, bed_name)
                print(f"Finished annotating {i}!\n")
            if args.classify:
                sat_db = find_files_with_suffix(args.classify,".bin")
                classify_seq(fasta_name,bed_name,args.kmer,sat_db)

    if args.mask:
        for fasta_file in seq_list:
            fasta_name = fasta_file[1]
            masked_fasta_name = f"{fasta_name}".split(".")[0] + "_masked.fa"
            print(masked_fasta_name)
            mask_fasta_with_bed(fasta_name, fasta_file[0], bed_name, masked_fasta_name)


if __name__ == "__main__":
    main()
