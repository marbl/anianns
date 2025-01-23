#!/usr/bin/env python3
import sys
import argparse
import math
import json
import os
import numpy as np
from moddotplot.parse_fasta import readKmersFromFile, getInputHeaders, isValidFasta
from moddotplot.estimate_identity import createSelfMatrix
from anianns.parse_input import readKmersFromFileNonHashed
from anianns.parse_dotplot import mask_fasta_with_bed, row_non_zero_stats, save_to_file, find_non_zero_length, mask_offdiagonals, mask_diagonals, merge_coordinates, create_bed_file
from collections import Counter

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
        required=True
    )
    parser.add_argument(
        "-b",
        "--band",
        type=float,
        default=8.0,
        help="Max height in Mbp of band."
    )
    parser.add_argument(
        "-k",
        "--kmer",
        default=21,
        help="k-mer length"
    )
    parser.add_argument(
        "--overlap",
        type=float,
        default=0.1,
        help="Percent overlap. Must be < 0.5."
    )
    parser.add_argument(
        "--identity",
        type=int,
        default=86,
        help="Identity threshold."
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default=None,
        help="Directory name for saving matrices and coordinate logs. Defaults to working directory.",
    )
    parser.add_argument(
        "-w",
        "--window",
        type=int,
        default=2000,
        help="Window size of ModDotPlot."
    )
    parser.add_argument(
        "-m",
        "--mask",
        action="store_true",
        help="Create a masked fasta file."
    )
    return parser

def main():
    args = get_parser().parse_args()
    ASCII_ART = """
     _          _      _                    
    / \   _ __ (_)    / \   _ __  _ __  ' ___ 
   / _ \ | '_ \| |   / _ \ | '_ \| '_ \  / __|
  / ___ \| | | | |  / ___ \| | | | | | | \__ |
 /_/   \_\_| |_|_| /_/   \_\_| |_|_| |_| |___/
    """
    print(ASCII_ART)
    # -----------INPUT SEQUENCE VALIDATION-----------
    seq_list = []
    for i in args.fasta:
        isValidFasta(i)
        headers = getInputHeaders(i)
        if len(headers) > 1:
            print(f"File {i} contains multiple fasta entries. \n")
            counter = 1
            for j in headers:
                counter += 1
        for j in headers:
            seq_list.append(j)
    if len(seq_list) > 1:
        print(f"Annotating the following sequences: {seq_list}. \n")
    else:
        print(f"Annotating {seq_list[0]}: \n")

    # -----------CREATE PATH IF NOT GIVEN------------
    if not args.output_dir:
        args.output_dir = os.getcwd()
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # -----------Run Ani Ann's on sequences-----------
    bed_coordinates = []
    for i in args.fasta:
        # Initialize bed file db
        k = readKmersFromFile(i, args.kmer, False)
        genome_length = len(k[0])
        divisor = math.ceil(genome_length/(args.band * 1000000))
        for j in range(0,divisor):
            start = int(j*1000000*args.band)
            end = int(min(start + (1000000 * args.band), genome_length))
            print(f"Creating a sub-matrix ({j+1}/{divisor}) of {seq_list[0]} from {start} to {end}: \n")
            submatrix = createSelfMatrix(len(k[0][start:end]),k[0][start:end],args.window,8,0.5,21,args.identity,False,500)
            nonzerodata = find_non_zero_length(submatrix,args.window)

            tuple_counts = Counter(nonzerodata[1])
            merged_result = merge_coordinates(tuple_counts)
            for (x, y), count in merged_result:
                x_coord = int((j*1000000*args.band) + x)
                y_coord = int((j*1000000*args.band) + y)
                #print(f"({x_coord}, {y_coord}): {count}")
                bed_coordinates.append((seq_list[0],x_coord,y_coord))
    create_bed_file(bed_coordinates, "anianns.bed")

    if args.mask:
        for i in args.fasta:
            masked_fasta_name = f"{i}".split(".")[0] + "_masked.fa"
            mask_fasta_with_bed(i,"tmp.bed",masked_fasta_name)

if __name__ == "__main__":
    main()
