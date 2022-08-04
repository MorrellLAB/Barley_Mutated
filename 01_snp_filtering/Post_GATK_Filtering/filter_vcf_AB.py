#!/usr/bin/env python3
"""This script filters by allele balance calculated from the GATK AD format field. This script sets samples to missing based on the thresholds. Afer running this script,
filter on proportion missing with bcftools. This script writes the filtered VCF lines to
standard output.

Usage: ./filter_vcf_AB.py [vcf_fp] [min_dev]

Where:
1. [vcf_fp] if the full filepath to the VCF file. [integer]
2. [min_dev] is the min amount of deviation from the expected 50-50 ref/alt reads. [proportion between 0 and 1]

This script was adapted from Li's script:
https://github.com/MorrellLAB/Barley_Inversions/blob/master/analyses/GATK_SNP_call/scripts/HeterozogotesVcfFilter.py.

This script can handle both phased and unphased genotypes.
"""

import sys
import os
import gzip

# User defined input arguments
# VCF filepath
vcf_fp = os.path.expanduser(sys.argv[1])
# Minimum amount of deviation from the expected 50-50 ref/alt reads.
# Example: 0.10 threshold would mean +/- 10% deviation from 50-50
min_dev = float(sys.argv[2])


def filter_vcf(f, min_dev):
    for line in f:
        #   Skip the header lines - write them out without modification
        if line.startswith('#'):
            sys.stdout.write(line)
        else:
            tmp = line.strip().split('\t')
            #   we aren't confident in our ability to call ancestral state of
            #   indels
            if len(tmp[3]) != 1 or len(tmp[4]) != 1:
                continue
            else:
                # VCFv4.2 format field: GT:AD:DP:GQ:PL
                format = tmp[8]
                sample_information = tmp[9:]
                #   The genotype information is the first element of each sample
                #   info block in the VCF
                for i, elem in enumerate(sample_information):
                    #   For the GATK HaplotypeCaller, the per-sample information
                    #   follows the form generally, but sometimes it has GAT:AD:DP:GQ:PGT:PID:PL, and majority of the time are only
                    #   GT:AD:DP:GQ:PL
                    info = elem.split(':')
                    gt = info[0] # Genotype information
                    #   We have to check for missing data first, because if it is
                    #   missing, then the other fields are not filled in
                    if not ('.' in gt or './.' in gt or '.|.' in gt):
                        # If genotype isn't missing, the AD field should be filled in
                        ad = info[1].split(',')
                        # GT is heterozygote, calculate allele balance and set to missing
                        #   if deviation is greater than threshold
                        if gt == '0/1' or gt == '0|1' or gt == '1|0':
                            ref = float(ad[0])
                            alt = float(ad[1])
                            if ref + alt != 0:
                                # AB defined the same way as Pedersen et al. 2021
                                balance = alt/(alt+ref)
                                # If the deviation is more than the threshold
                                # set the genotype to missing
                                if abs(0.5 - balance) > min_dev:
                                    info[0] = './.'
                                    reformated_info = ':'.join(info)
                                    sample_information[i] = reformated_info
                # Append new sample information to first 9 fields
                updated_sample_info = tmp[:9] + sample_information
                sys.stdout.write('\t'.join(updated_sample_info) + '\n')

#   Read the file in line-by-line
if "gz" in vcf_fp:
    with gzip.open(vcf_fp, 'rt') as f:
        filter_vcf(f, min_dev)
else:
    with open(vcf_fp, 'rt') as f:
        filter_vcf(f, min_dev)
