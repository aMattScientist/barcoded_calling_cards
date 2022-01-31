This is a fork of Arnav Moudgil's Calling Card analysis repository (https://github.com/arnavm/calling_cards), with minor modifications to enable analysis of barcoded self-reporting transposon calling cards. Reference: https://www.biorxiv.org/content/10.1101/2021.04.15.439516

# Prerequisites

The scripts in this folder require `python3`, preferrably a relatively recent version (e.g. ≥ 3.5). In addition, you will need to install the following modules:
- `pysam`
- `numpy`
- `pandas`
- `twobitreader`
- `pybedtools`
- `scipy`
- `statsmodels`
- `astropy`

The easiest way to install them is:
- If `python3` is the default:
`pip install pysam numpy pandas twobitreader pybedtools scipy statsmodels astropy`
- Otherwise:
`pip3 install pysam numpy pandas twobitreader pybedtools scipy statsmodels astropy`

To figure out which version of `python` is installed by default:
- `python -V`
