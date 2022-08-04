#!/usr/bin/env python3
"""A script to convert a VCF to tab-delimited output compatible with
bcftools annotate.

This script adds the annotations:
1) novel
2) known

Note: novel vs known variants here are solely based on comparing to a vcf of known variants

Usage: ./vcf_to_bcftools_ann_tsv-known.py [vcf] [known_variants_vcf] > out_ann.txt

Where:
1) [vcf] is the vcf file containing variants we want to annotate
2) [known_variants_vcf] is the vcf containing a set of known variants
    (e.g., for barley those could be BOPA, 9K, or 50K SNPs)
"""

import sys
from cyvcf2 import VCF

# Print usage message if we are missing input arguments
if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(0)

# User provided command line arguments
vcf_fp = sys.argv[1]
known_vcf_fp = sys.argv[2]

# Load known variants vcf
# Store in dictionary for fast lookup
known_dict = {}
for variant in VCF(known_vcf_fp):
    record = '_'.join([variant.CHROM, str(variant.POS)])
    # Only need to lookup the keys, so we don't need the other values here
    known_dict[record] = []

# Load vcf and prepare column values
# bcftools annotate requires the following columns:
#   CHROM, POS, REF, ALT, QUAL, and NEW_TAG
for variant in VCF(vcf_fp):
    record = '_'.join([variant.CHROM, str(variant.POS)])
    if record in known_dict.keys():
        # This is a known variant, add appropriate annotation
        print('\t'.join([variant.CHROM, str(variant.POS), variant.REF, ','.join(variant.ALT), str(variant.QUAL), "known"]))
    else:
        # This variant doesn't exist in the known variants VCF, so we'll
        #   consider it a novel variant
        print('\t'.join([variant.CHROM, str(variant.POS), variant.REF, ','.join(variant.ALT), str(variant.QUAL), "novel"]))
