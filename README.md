![Release](https://img.shields.io/github/v/release/marbl/anianns?sort=semver&label=stable%20release)
![Coverage](https://img.shields.io/codecov/c/github/marbl/anianns?label=coverage)
![Black](https://img.shields.io/badge/code%20style-black-000000.svg)

![](images/anianns_logo.png)

- [About](#about)
- [Installation](#installation)
- [Usage](#usage)
  - [Annotate Mode](#annotate-mode)
    - [Required args](#required-args)
    - [Optional args](#optional-args)
    - [Sample run](#sample-run)
    - [Repeat Masking](#repeat-masking)
  - [Build db](#build-db)
    - [Required args](#required-args-1)
    - [Optional Arguments](#optional-arguments)
    - [Sample database](#sample-database)
  - [NTRPrism](#ntrprism)
- [Questions](#questions)


## About
_AniAnn's_: __ANI__ Inferred **ANN**otation of Tandem Repeats

_AniAnn's_ is an _a priori_ satellite array detection and annotation software package. _AniAnn's_ uses a matrix of Average Nucleotide Identity (ANI) values similar to [ModDotPlot](https://github.com/marbl/ModDotPlot) in order to infer the location and orientation of satellite arrays. It introduces new downstream analysis to accurately annotate its contents and boundaries.

## Installation

You can download the current release from GitHub by using: 

```
git clone https://github.com/marbl/anianns.git
cd anianns
```

Although optional, we recommend setting up a virtual environment:

```
python -m venv venv
source venv/bin/activate
```

Once the virtual environment is activated, you can install the required dependencies:

```
python -m pip install .
```

By default, _AniAnn's_ installs **without** _ModDotPlot_ as a dependency. If you would like to include _ModDotPlot_ into the same venv for plotting, you can install using the following command:

```
python -m pip install .[moddotplot]
```

_AniAnn's_ is also available to install using PyPI:

```
pip install anianns
```

Once installed, confirm _AniAnn's_ was installed correctly by running `python -m anianns -h`, or simply with the shortcut `anianns -h`:

```

                                                 __       __
     _          _     _                        .'  `'._.'`  '.
    / \   _ __ (_)   / \   _ __  _ __  ' ___  |  .--;   ;--.  | 
   / _ \ | '_ \| |  / _ \ | '_ \| '_ \  / __| |  (  /   \  )  | 
  / ___ \| | | | | / ___ \| | | | | | | \__ |  \  ;` /^\ `;  / 
 /_/   \_\_| |_|_|/_/   \_\_| |_|_| |_| |___/   :` .'._.'. `;
                                                '-`'.___.'`-'

usage: anianns [-h] {annotate,build_db} ...

Ani Ann's: ANI Inferred ANNotation of Tandem Repeats

positional arguments:
  {annotate,build_db,ntrprism}  Choose mode: annotate, build_db, ntrprism
    annotate           Takes an input fasta and outputs an annotated bedfile of satellite arrays.
    build_db           Takes an input fasta, a bedfile of satellite coordinates, and an optional config file and outputs a db of satellite k-mers.
    ntrprism           Takes an input fasta, and a bedfile or region. Outputs k-mer histogram.


options:
  -h, --help       show this help message and exit
```

Note that _AniAnn's_ might take a while to run during your first usage. This is because the Python interpreter is compiling source code into the pycache directory. Subsequent runs will use the pre-compiled code and load much faster!

## Usage 

_AniAnn's_ must be run either in `annotate` mode, or `build_db` mode. 

### Annotate Mode

```
anianns annotate -f <FASTA_FILENAME(S)> <ARGS>
```

_AniAnn's_ requires at least one FASTA file as input. It generates one output annotation file (default BED) per sequence contained in the input FASTA file(s). Output annotation files are named based on the sequence identifier from the FASTA header.

Annotation of arrays into known satellite classes **must** be done through the use of a satellite _k_-mer database, using the command `--classify <directory>`. See [creating an annotation database](#creating-an-annotation-database) for more information.

#### Required args

`-f / --fasta <FILENAME(S)>`

Fasta file(s) to input. Multifasta files are accepted. 

#### Optional args

`-s / --seq_id <STR>`

Sequence ID to extract (multiple if using multifasta file). Will ignore if not found. **Default: None.**

`-d / --directory <DIR>`

Name of output directory. **Default: current working directory.**

`-o / --output-format <STR>`

Output annotation file format. Options are BED, GTF, GFF, CSV, TSV, JSON. **Default: BED.**

`-m / --mask <ALL>`

Name or repeat length of satellite arrays to mask. Replaces deteced satellites with N's. See [Repeat Masking](#repeat-masking) for more info. **Default: None.**

`--soft <BOOL>`

Softmask flag. Instead of N's, will force bases lower-case in detected satellite arrays. Must be used with `--mask`. **Default: None.**

`-c / --classify <DIR>`

Directory containing .db or .msh k-mer db files. Required for annotation into known satellites. **Default: None.**

`-t / --threshold <INT>`

Confidence threshold. This relates to the minimum percentage of _k_-mers in a satellite array required to be within that are contained within the selected. A lower number will be more sensitive, at the risk  **Default: 50.**

`-k / --kmer <INT>`

K-mer size to use. This should be large enough to distinguish unique k-mers with enough specificity, but not too large that sensitivity is removed. **Default: 21.**

`-i / --identity <INT>`

Minimum sequence identity cutoff threshold when running ModDotPlot. While it is possible to go as low as 50% sequence identity, anything below 80% is not recommended. **Default: 86.** 

`-w / --window <INT>`

The central dotplot window size. AniAnn's automatically scans half, the supplied
value, and double the supplied value. For example, `-w 5000` scans windows of
2500, 5000, and 10000 bp. AniAnn's hashes each sequence band once, scans the
additional resolutions without dense matrices, reconciles overlapping calls,
and boundary-refines only the winning call at each locus. It also writes a
`<sequence>_window_selection.tsv` audit file recording the chosen resolution.
Comparable spans favor the finer resolution; a coarser call wins when it
recovers materially more supported array sequence.
**Default: 2000 (scans 1000, 2000, and 4000).**

`--band <FLOAT>`

_AniAnn's_ streams the sequence in bands of this many megabases. Normal
annotation scans outward from the diagonal without materializing a full band
matrix; `--plot` or `--distal` creates the primary-resolution matrix needed for
visualization or distal-link detection. **Default: 2.**

`--cache-dir <DIR>`

Directory for reusable canonical k-mer hashes. By default, AniAnn's uses a
persistent per-user cache shared by every output directory. Set
`ANIANNS_CACHE_DIR` to change that shared location globally, or use this option
for one run.

Every annotation run also writes `satellite_dsu.tsv` and a human-readable
`satellite_dsu.txt`. Each boundary-refined
satellite is a DSU node, including unlinked singletons. With `--distal`,
validated distal links
union their two satellite nodes, so chains of distal relationships share one
stable `component_id`; the table also reports component size and direct-link
count.

`--identifier <STR>`

Name of identifier. Used when no matches to a k-mer db are found, or if `--classify` is not provided. bed file to output to. **Default: None.**

`-p / --plot <bool>`

Save a low-resolution PNG heatmap for every `--band` length matrix. Heatmaps are
collected under `<output directory>/matrix_plots` and are not displayed
interactively. Plotting alone does not search for distal links and does not add
red overlays. **Default: disabled.**

`--distal`

Detect distal-satellite links independently of whether plots are requested.
Detected off-diagonal blocks are written to
`<sequence>_distal_links.bedpe`, and used to link satellites in the DSU. A
compact cross-band matrix and its genomic coordinates are saved to
`<sequence>_distal_neighborhood.npz`; two neighboring windows on either side of
each candidate are retained. An unmatched distal axis must pass NTR Prism and
boundary refinement before it is added to BED/CSV or linked in the DSU.
Neighboring bands are bridged with bounded candidate-to-all comparisons in both
directions plus a 50-window seam scan. With both `--distal --plot`, links are
outlined in red. When the bridge finds evidence, a sparse
two-band heatmap is saved under `matrix_pairs/`; adjacent pairs are then skipped
by the final global comparison to avoid duplicate work.
**Default: disabled.**

`--verbose <bool>`

Verbose logging output. Creates a log file at `--directory`. **Default: None.**

`--quiet <bool>`

Suppress all logging output. **Default: None.**

#### Sample run

`anianns -f sample_sequences/sample_hap1_.fa`

Upon running the above command, you should see the following output in `sample_hap1.bed`:

```
#chrom	start	end	name	score	strand	thickStart	thickEnd	itemRgb
sample_hap1	67078	407613	147	.	67078	407613	230,57,70
sample_hap1	634567	669509	5	.	634567	669509	42,157,143
sample_hap1	669978	1148691	68	.	669978	1148691	241,196,15
sample_hap1	1148878	1472729	5	.	1148878	1472729	42,157,143
sample_hap1	1472778	2702615	42	.	1472778	2702615	138,43,226
sample_hap1	2703054	2729531	68	.	2703054	2729531	241,196,15
sample_hap1	2729678	2834698	5	.	2729678	2834698	42,157,143
sample_hap1	3044365	3208742	147	.	3044365	3208742	230,57,70
sample_hap1	3225177	3297021	48	.	3225177	3297021	30,144,255
```

This is a BED file containing inferred satellite intervals. The value in the `score` column indicates the periodicity of the satellite array. Since no k-mer database was used, _AniAnn's_ does not attempt to classify each array. _AniAnn's_ will label arrays it determines to be related as the same color in `itemRgb`.

To classify each line of the BED file into a known satellite array, a database of _k_-mers must be used. See [creating an annotation database](#creating-an-annotation-database) for more information. We will use the following provided _k_-mer db for our sample run:

`anianns -f sample_sequences/sample_hap1.fa --classify sample_sequences/sample_kmer_db`

```
#chrom	start	end	name	score	strand	thickStart	thickEnd	itemRgb
sample_hap1	67078	407613	ACRO	147	.	67078	407613	230,57,70
sample_hap1	634567	669509	HSat3	5	.	634567	669509	42,157,143
sample_hap1	669978	1148691	bSat	68	.	669978	1148691	241,196,15
sample_hap1	1148878	1472729	HSat3	5	.	1148878	1472729	42,157,143
sample_hap1	1472778	2702615	HSat1A	42	.	1472778	2702615	138,43,226
sample_hap1	2703054	2729531	bSat	68	.	2703054	2729531	241,196,15
sample_hap1	2729678	2834698	HSat3	5	.	2729678	2834698	42,157,143
sample_hap1	3044365	3208742	ACRO	147	.	3044365	3208742	230,57,70
sample_hap1	3225177	3297021	CER	48	.	3225177	3297021	30,144,255
```

Note that you *must* use the default k-mer value as the classification database. The default _k_ = 21 

#### Repeat Masking

```
anianns annotate -f <FASTA_FILENAME(S)> --mask <ARGS>
```

Use of _AniAnn's_ as a tool to mask satellite arrays from a given sequence is done through the `-m/--mask` argument. _AniAnn's_ will output each sequence into its own masked fasta file in the `--directory` ouptut folder. By default, running `--mask` with no parameters will mask everything deemed a satellite array.

- If string parameter(s) are provided (e.g., `hSat1 hSat2`), _AniAnn's_ will mask arrays that match the provided class. Note this **must** be used in conjunction with `--classify` in order to match names. Matches are case-insensitive.
- If integer parameter(s) are provided (e.g., `6 7 42`), _AniAnn's_ will mask arrays whose predominant monomer length is equal to that provided. 


### Build db

```
anianns build_db -f <FASTA_FILENAME(S)> -b <BED_FILENAME(S)> <ARGS>
```

Building a database of satellite _k_-mers is how _AniAnn's_ can match a detected satellite array into a known repeat class. _AniAnn's_ will extract and compress _k_-mers at the coordinates provided by a bed file. Running this will create a directory at `--directory`. 

#### Required args

`-f / --fasta <FILENAME(S)>`

Fasta file(s) to input. Multifasta files are accepted. 

`-b / --bed <FILENAME(S)>`

Bed file(s) to input. Must contain at least 4 columns: 1 chrom, 2 start, 3 end, 4 name.  

#### Optional Arguments

`-c / --config <FILENAME>`

Name of config file to use. See [Sample database for more info](#sample-database)**Default: None.**

`-k / --kmer <INT>`

K-mer size. Note that k-mers **must** be the same length when running `anianns annotate` in order for classification to work. **Default: 21.**

`-d / --directory <DIR>`

Specifies the output directory where results will be saved. If not provided, _AniAnn's_ automatically creates an output folder in the current working directory:

- If the BED file contains a track header, the folder name will be derived from that header.
- If no header is present, the folder name will default to the k-mer length (e.g., k_21).
**Default: current working directory.**

`-v / --verbose <BOOL>`

Verbose logging output. Will output a log file into `--directory`. **Default: False.**

`-q / --quiet <BOOL>`

Suppress all logging output. **Default: False.**

#### Sample database

Here, we will create a db of known satellite arrays using the HG002 human genome CenSat annotation track ([credit: Hailey Loucks](https://github.com/hloucks/CenSatData/tree/main)):

`wget https://raw.githubusercontent.com/hloucks/CenSatData/refs/heads/main/HG002/v1.1/hg002v1.1.cenSatv2.0.bed`

We want to remove any centromere transition regions `ct` from this BED file:

`awk 'NR==1 || $4!="ct"' hg002v1.1.cenSatv2.0.bed > hg002v1.1.cenSatv2.0.ctRemoved.bed`

Finally, we want to group related satellite arrays into the same class. This is done using a `--config` file. For this example, we will use the file provided in `config` folder of this repo: 

```
head config/sample_config.json 
{
    "gSat": [
      "gSat(GSAT,GSATX)",
      "gSat(GSAT)",
      "gSat(GSATII,GSATX)",
      "gSat(GSATII,TAR1)",
      "gSat(GSATII)",
      "gSat(GSATX)",
      "gSat(TAR1)"
    ],
```

This merges the _k_-mers of all variations of gSat into the same class. Anything not in the config file is not included in the *k*-mer db. The config file must be in standard JSON format. If a config file is not provided, each unique value in column 4 of the input bed file will become its own unique class. 

Creating a *k*-mer db for HG002 using the provided config file takes around 3 minutes. This results in 16.9 million unique _k_-mers, compressed down into a 53mb directory. Note that increasing the *k*-mer size will increase the directory size, as a more specific *k*-mer threshold will increase the total number of unique *k*-mers.

### NTRPrism

NTRPrism reports the most common distances between consecutive occurrences of
the same forward k-mer in one FASTA region. Nearby spacing values are combined
before ranking and plotting, so values such as 170 and 171 bp contribute to the
same peak by default.

```bash
anianns ntrprism \
  -f assembly.fa.gz \
  -s chr1_MATERNAL \
  -r 120000000 121000000 \
  -k 6
```

The region uses 0-based, half-open coordinates. By default, the command prints
the ten strongest merged spacing peaks, each peak's count as a percentage of
the requested interval, and a horizontal ASCII histogram. It does not create
files unless `--save` is supplied:

```bash
anianns ntrprism \
  -f assembly.fa.gz \
  -s chr1_MATERNAL \
  -r 120000000 121000000 \
  -k 11 \
  --save \
  -d chr1_ntrprism
```

With `--save`, AniAnn's also writes a top-ten text report and a PNG histogram.
The `-k / --kmer` option controls the k-mer length and defaults to 21, matching
the annotation pipeline. The
default nearby-value merge distance is 1 bp; use `--merge-distance` to change
it. Use `--quiet --save` for files without the terminal report. Missing or
unreadable FASTA files, unknown sequence IDs, invalid ranges, and regions
shorter than the selected k-mer size produce an error and a nonzero exit status.

## Questions

For bug reports or general usage questions, please raise a GitHub issue, or email alex ~dot~ sweeten ~at~ nih ~dot~ gov
