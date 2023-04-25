#!/usr/bin/env python3
"""Per sample counts for indels in VCF.
Also returns columns for per sample counts: total per sample snps, total het, total hom_alt.

Get counts for the following:
1) num_snps: number of het and hom_alt SNPs per individual
2) num_hom_alt
3) num_het
4) prop_hom_alt
5) prop_het

# Usage: ./count_indels_per_individual.py [vcf_fp] [out_dir] [out_prefix]
"""

import sys
import os
import pandas as pd
import numpy as np
from cyvcf2 import VCF

if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(0)

# User provided input arguments
# VCF containing all samples
vcf_fp = sys.argv[1]
# Output directory
out_dir = os.path.expanduser(sys.argv[2])
# Output file prefix
out_prefix = sys.argv[3]

# Note: not very pretty but works and gets the info needed!
# Read VCF
vcf = VCF(vcf_fp)
vcf_list = [['id', 'chr', 'pos'] + vcf.samples]

# Label genotypes as miss, hom_ref, het, and hom_alt
# Note: current cyvcf2 gt_types function doesn't handle missing genotypes correctly
# gt_types incorrectly categorizes ./. genotypes as HOM_ALT, hence, the alternative solution below
for variant in vcf:
    gt_list = []
    for gt in variant.genotypes:
        if gt[0] == -1 and gt[1] == -1:
            gt_list.append('miss')
        elif gt[0] != -1 and gt[1] != -1:
            if gt[0] == 0 and gt[1] == 0:
                gt_list.append('hom_ref')
            elif gt[0] != gt[1]:
                gt_list.append('het')
            elif gt[0] != 0 and gt[1] != 0 and gt[0] == gt[1]:
                gt_list.append('hom_alt')
    vcf_list.append(['_'.join([variant.CHROM, str(variant.POS)]), variant.CHROM, variant.POS] + gt_list)
# Convert to pandas dataframe and set first row as header
vcf_df = pd.DataFrame(vcf_list)
vcf_df = vcf_df.T.set_index(0).T

# Per sample counts
counts_summary = [['sample', 'num_indels', 'num_hom_alt_indels', 'num_het_indels', 'prop_hom_alt_indels', 'prop_het_indels']]
for column in vcf_df:
    if column != "id" and column != "chr" and column != "pos":
        counts = vcf_df[column].value_counts()
        num_snps = counts['hom_alt']+counts['het']
        counts_summary.append([column, num_snps, counts['hom_alt'], counts['het'], counts['hom_alt']/num_snps, counts['het']/num_snps])

# Convert to pandas dataframe and set first row as column names
counts_summary_df = pd.DataFrame(counts_summary)
counts_summary_df = counts_summary_df.T.set_index(0).T

# Save final counts to table
counts_summary_df.to_csv(out_dir + '/' + out_prefix + '_per_sample_counts.txt', sep='\t', header=True, na_rep='NA', index=False)
