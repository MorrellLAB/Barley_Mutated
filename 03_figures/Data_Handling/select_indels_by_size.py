#!/usr/bin/env python3
"""Filter indels by minimum and maximum sizes inclusive.

Usage: ./select_indels_by_size.py [vcf_fp] [min_size_bp] [max_size_bp] > output.vcf
"""

import sys
import os
from cyvcf2 import VCF

if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(0)

# User provided input argument
vcf_fp = os.path.expanduser(sys.argv[1])
min_size = int(sys.argv[2])
max_size = int(sys.argv[3])

def read_indels(vcf_fp, min_size, max_size):
    vcf = VCF(vcf_fp)
    # Print header to stdout
    print(str(vcf.raw_header).strip())
    for v in vcf:
        # Insertions will be a positive number and deletions will be negative
        # Take absolute value for easier filtering
        indel_size = abs(len(v.ALT[0])-len(v.REF))
        if indel_size >= min_size and indel_size <= max_size:
            # Convert back to string and print
            print(str(v).strip())


read_indels(vcf_fp, min_size, max_size)
