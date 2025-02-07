#!/usr/bin/env python
# coding: utf-8

# This program takes as input an aligned BAM file of calling card reads 
# It then tags each read with the SRT barcode that was appended to read name by umi-tools 
# Prefix for SRT barcode = XQ

import argparse
import pysam
import sys

parser = argparse.ArgumentParser()
parser.add_argument("input", type = str)
parser.add_argument("output", type = str)
parser.add_argument("-t", "--tag", type = str, required = True)
args = parser.parse_args()

if __name__ == "__main__":
    # Open files
    bamfile = pysam.AlignmentFile(args.input, "rb")
    outfile = pysam.AlignmentFile(args.output, "wb", template = bamfile)
    tag = args.tag.split(':')
    for read in bamfile:
        srt_bc = "".join(read.query_name.split('_')[-1:]) ## Last element of the read name is the SRT barcode from UMI tools (4 bp CellBarcode and 0 bp UMI)
        read.set_tag(tag[0], srt_bc, tag[1])
        outfile.write(read)
    # Close files
    outfile.close()
    bamfile.close()
