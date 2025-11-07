import struct
ASCII_ART = """
                                                 __       __
     _          _     _                        .'  `'._.'`  '.
    / \   _ __ (_)   / \   _ __  _ __  ' ___  |  .--;   ;--.  | 
   / _ \ | '_ \| |  / _ \ | '_ \| '_ \  / __| |  (  /   \  )  | 
  / ___ \| | | | | / ___ \| | | | | | | \__ |  \  ;` /^\ `;  / 
 /_/   \_\_| |_|_|/_/   \_\_| |_|_| |_| |___/   :` .'._.'. `;
                                                '-`'.___.'`-'
"""
BED_COLUMNS = [
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "itemRgb",
    "blockCount",
    "blockSizes",
    "blockStarts",
]
DESCRIPTION = "Ani Ann's: ANI Inferred ANNotation of Tandem Repeats"
HEADER_FORMAT = "32sQ"  # 32-byte name (padded), 8-byte unsigned int (number of kmers)
HEADER_SIZE = struct.calcsize(HEADER_FORMAT)
VERSION = "0.5.0"

