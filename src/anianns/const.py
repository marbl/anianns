import struct
VERSION = "0.3.0"
HEADER_FORMAT = "32sQ"  # 32-byte name (padded), 8-byte unsigned int (number of kmers)
HEADER_SIZE = struct.calcsize(HEADER_FORMAT)

COLOR_MAPPING = {
    'hsat1a':   '0,222,96',
    'hsat1b':   '200,222,96',
    'hsat2':    '31,81,245',
    'hsat3':    '51,81,137',
    'sst1': '172,51,199',
    'dhor': '255,146,0',
    'mixedalpha': '255,146,0',
    'hor': '255,146,0',
    'active_hor': '153,0,0',
    'bsat': '250,153,255',
    'gsat': '180,153,255',
    'other':    '128,128,128'
}