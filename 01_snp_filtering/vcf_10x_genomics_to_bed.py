#!/usr/bin/env python3
"""Script that converts 10x Genomics formatted VCF file to BED format.
This script handles non-BND variants that have end positions stored
in the INFO/END field.

Usage: ./vcf_10x_genomics_to_bed.py [vcf_file] > out_file.bed

Where:
1) [vcf_file] is the full filepath to the VCF file
"""

import sys
import os

if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(0)

# User provided input arguments
VCF = os.path.expanduser(sys.argv[1])

with open(VCF, 'r') as f:
    for line in f:
        if line.startswith("#"):
            continue
        else:
            vl = line.strip().split('\t')
            chrom = vl[0]
            # BED is 0 indexed
            start = str(int(vl[1])-1)
            end = vl[7].split(';')[0].split('=')[1]
            print('\t'.join([chrom, start, end]))
