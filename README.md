This is a fork of Arnav Moudgil's Calling Card analysis repository (https://github.com/arnavm/calling_cards), with minor modifications to enable analysis of barcoded self-reporting transposon calling cards. Reference: https://www.biorxiv.org/content/10.1101/2021.04.15.439516

# Prerequisites

The scripts in this folder require `python3`, preferrably a relatively recent version (e.g. ≥ 3.5). In addition, you will need to install the following modules:
- `pysam`
- `numpy`
- `pandas`
- `twobitreader`
- `pybedtools`
- `scipy`

The easiest way to install them is:
- If `python3` is the default:
`pip install pysam numpy pandas twobitreader pybedtools scipy`
- Otherwise:
`pip3 install pysam numpy pandas twobitreader pybedtools scipy`

To figure out which version of `python` is installed by default:
- `python -V`

Running this software also requires several common bioinformatics modules, ideally on a Linux/Unix based high-performance computing environment. 

The main shell script requires: 
- `cutadapt`
- `umi-tools`
- `star or (NovoAlign)` 
- `samtools`
- `bedtools`
- `java`

# Getting Started

Create a new directory, for example, named CALLINGCARDS.  
Within this directory, make 3 more directories:  output_and_analysis, raw, and CODE.  
Download the contents of this repository to the CODE directory.  





