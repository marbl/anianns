![](images/anianns_logo.png)

- [About](#about)
- [Installation](#installation)
- [Usage](#usage)
- [Questions](#questions)


## About

Ani Ann's is an _a priori_ satellite detection and annotation software. Ani Ann's uses a matrix of Average Nucleotide Identity (Ani) values created by [ModDotPlot](https://github.com/marbl/ModDotPlot) to infer the location and type of satellites. 

Note that Ani Ann's is currently under active development, and not all features are currently available.

## Installation

```
git clone https://github.com/marbl/ModDotPlot.git
cd ModDotPlot
```

Although optional, setting up a virtual environment is recommended:

```
python -m venv venv
source venv/bin/activate
```

Once activated, you can install the required dependencies:

```
python -m pip install .
```

## Usage

Currently, use of is limited to . Classification of satellites into . 

`-f / --fasta <file>`

Fasta files to input. Multifasta files are accepted. 

`-k / --kmer <int>`

K-mer size to use. This should be large enough to distinguish unique k-mers with enough specificity, but not too large that sensitivity is removed. Default: 21.

`-o / --output-dir <string>`

Name of output directory. Default is current working directory.

`-id / --identity <int>`

Minimum sequence identity cutoff threshold when running ModDotPlot. Default is 86. While it is possible to go as low as 50% sequence identity, anything below 80% is not recommended. 

`-w / --window <int>`

Dotplot window size, or the number of bp contained within each pixel in a plot. This is proportional to the sensitivity of satellite detection (ie. lower is more accurate, at the expense of runtime). Default is 2000.

`-b / --band <int * 1000000>`

When creating dotplots, to save time, multiple plots of a certain size are used instead of the entire seqeunce length. This can be adjusted here (default: 8, or 8Mbp). Increasing this will improve detection of off-target satellites, at the expense of runtime.

`--bed <str>`

Name of bed file to output to. Default is `anianns-output.bed`. 

## Questions

For bug reports or general usage questions, please raise a GitHub issue, or email alex ~dot~ sweeten ~at~ nih ~dot~ gov
