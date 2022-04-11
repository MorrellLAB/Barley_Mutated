#!/usr/bin/env python3
"""Script that converts PacBio/Nanopore VCF called with Sniffles2 to BED format.
This script handles non-BND variants that have end positions stored
in the INFO/END field. BND variants are skipped (usually on different chromosomes
and for our purposes can be hard to visualize).

Usage: ./vcf_sniffles_long_read_to_bed.py [vcf_file] > out_file_noBND.bed

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
            # Only converting to BED for non-BND variants
            # Skip BND variant types
            if 'BND' not in vl[2]:
                chrom = vl[0]
                # BED is 0 indexed
                start = str(int(vl[1])-1)
                info_field = vl[7].split(';')
                end_pos_idx = [i for i, elem in enumerate(info_field) if 'END=' in elem]
                end = info_field[end_pos_idx[0]].split('=')[1]
                print('\t'.join([chrom, start, end]))
