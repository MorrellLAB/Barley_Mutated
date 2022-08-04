#!/usr/bin/env python3
"""A script to convert a VCF to tab-delimited output compatible with
bcftools annotate.

This script adds the annotations:
1) filtered - for variants that are filtered out
2) retained_filt1 - for variants that are retained after first pass filtering
3) retained_filt2 - for variants that are retained after second pass filtering

These annotations allow help us fine tune our filtering.

Usage: ./vcf_to_bcftools_ann_tsv-filtered.py [vcf] [filt1_vcf] [filt2_vcf] > out_ann.txt

Where:
1) [vcf] is the vcf file containing variants we want to annotate (e.g., unfiltered vcf)
2) [filt1_vcf] is the vcf containing variants retained after the first pass of filtering
3) [filt2_vcf] is the vcf containing variants retained after the second pass of filtering
"""

import sys
from cyvcf2 import VCF

# Print usage message if we are missing input arguments
if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(0)

# User provided command line arguments
vcf_fp = sys.argv[1]
filt1_vcf_fp = sys.argv[2]
filt2_vcf_fp = sys.argv[3]

# Load filt1 variants vcf
# Store in dictionary for fast lookup
filt1_dict = {}
for variant in VCF(filt1_vcf_fp):
    record = '_'.join([variant.CHROM, str(variant.POS)])
    # Only need to lookup the keys, so we don't need the other values here
    filt1_dict[record] = []

# Load filt2 variants vcf
filt2_dict = {}
for variant in VCF(filt2_vcf_fp):
    record = '_'.join([variant.CHROM, str(variant.POS)])
    # Only need to lookup the keys, so we don't need the other values here
    filt2_dict[record] = []

# Load vcf and prepare column values
# bcftools annotate requires the following columns:
#   CHROM, POS, REF, ALT, QUAL, and NEW_TAG
for variant in VCF(vcf_fp):
    record = '_'.join([variant.CHROM, str(variant.POS)])
    if (record not in filt1_dict.keys()) and (record not in filt2_dict.keys()):
        # Variant is filtered out, add appropriate annotation
        print('\t'.join([variant.CHROM, str(variant.POS), variant.REF, ','.join(variant.ALT), str(variant.QUAL), "filtered"]))
    elif (record in filt1_dict.keys()) and (record not in filt2_dict.keys()):
        # retained_filt1 but removed after second pass filtering, add appropriate annotation
        print('\t'.join([variant.CHROM, str(variant.POS), variant.REF, ','.join(variant.ALT), str(variant.QUAL), "retained_filt1"]))
    elif (record in filt1_dict.keys()) and (record in filt2_dict.keys()):
        # retained_filt2 (passed all levels of filtering), add appropriate annotation
        #   consider it a novel variant
        print('\t'.join([variant.CHROM, str(variant.POS), variant.REF, ','.join(variant.ALT), str(variant.QUAL), "retained_filt2"]))
