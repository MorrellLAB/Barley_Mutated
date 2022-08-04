#!/usr/bin/env python3
"""This script adds the allele balance (AB) annotation to the format field for
each sample. AB is calculated from the AD annotation output by GATK.

Allele balance is calculated in the same way as defined in Pedersen et al. 2021:
alternate reads / (alternate reads + reference reads)."""

import sys
import os
import gzip
import vcfpy

# User defined input arguments
#vcf_fp = os.path.expanduser(sys.argv[1])
vcf_fp = os.path.expanduser("~/Downloads/temp_msi/toy.vcf")

#reader = vcfpy.Reader.from_path(vcf_fp)

for record in reader:
    print(record.get_format_field_info('AD'))

for record in vcfpy.Reader.from_path(vcf_fp):
    for c in record.calls:
        #print(c)
        gt = c.data.get('GT')
        ad = c.data.get('AD')
        #print(ad)
        ref = ad[0]
        alt = ad[1]
        if gt == "./." or gt == ".|.":
            # AB set to missing
            tmp_ab = "."
        elif len(ad) > 2:
            # AB set to missing
            # There is more than one alternate allele present
            print(ad)
        elif alt + ref != 0:
            # Calculate
            tmp_ab = alt/(alt+ref)
        elif alt == 0 and ref == 0:
            # AB is set to missing
            tmp_ab = "."
