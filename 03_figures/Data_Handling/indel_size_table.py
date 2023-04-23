#!/usr/bin/env python3
"""For each indel size class, generate the proportion of indels in that class.
Generates indel size table for plotting purposes. Output file columns include:

1) indel_size
2) proportion
3) dataset

Requires Python >= 3.9

Usage: ./indel_size_table.py [vcf_fp] [out_dir] [dat_name]

Where:
1) [vcf_fp] is the full filepath to VCF
2) [out_dir] is the full output directory path
3) [dat_name] is the name to fill the "dataset" column in the output
"""

import os
import sys
from cyvcf2 import VCF
import pandas as pd
import numpy as np

if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(0)

# User provided input arguments
vcf_fp = os.path.expanduser(sys.argv[1])
out_dir = os.path.expanduser(sys.argv[2])
dat_name = sys.argv[3]

def read_indels(vcf_fp, dataset_name):
    vcf = VCF(vcf_fp)
    indel_sizes_list = []
    for v in vcf:
        # Insertions will be a positive number and deletions will be negative
        indel_size = len(v.ALT[0])-len(v.REF)
        if indel_size >= 0:
            dat_label = "ins"
        elif indel_size < 0:
            dat_label = "del"
        indel_sizes_list.append([v.CHROM, v.POS, v.REF, v.ALT[0], indel_size, dat_label, dataset_name])
    # Convert to pandas dataframe
    indel_sizes_df = pd.DataFrame(indel_sizes_list,
                                  columns=['chr', 'pos', 'ref', 'alt', 'indel_size', 'var_type', 'dataset'])
    # Calculate proportion in each size group
    prop_df = indel_sizes_df.groupby('dataset')['indel_size'].value_counts(normalize=True).fillna(0).unstack(-1).transpose().reset_index()
    # Rename columns
    prop_df.rename(columns={dataset_name: "proportion"}, inplace=True)
    # Add dataset category
    prop_df['dataset'] = dataset_name
    return(prop_df)


df = read_indels(vcf_fp, dat_name)

# Prepare VCF basename
vcf_bn = os.path.basename(vcf_fp)
# .removesuffix requires python >= 3.9
if '.vcf.gz' in vcf_bn:
    vcf_prefix = vcf_bn.removesuffix('.vcf.gz')
else:
    vcf_prefix = vcf_bn.removesuffix('.vcf')
# Prepare output file name from vcf file prefix
out_fp = out_dir + '/indel_sizes-' + vcf_prefix + '.txt'

# Save to file
# Exclude row names, keep header
df.to_csv(out_fp, sep='\t', index=False)
