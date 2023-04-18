#!/usr/bin/env python3
"""Per sample counts for "Deleterious" vs "Tolerated" and synonymous vs nonsynonymous.
Also returns columns for per sample counts: total per sample snps, total het, total hom_alt.

Get counts for the following:
1) num_snps: number of het and hom_alt SNPs per individual
2) num_hom_alt
3) num_het
4) prop_hom_alt
5) prop_het
6) het-Deleterious
7) het-Tolerated
8) hom_alt-Deleterious
9) hom_alt-Tolerated
10) total_del
11) total_tol
12) het-nonsyn
13) hom_alt-nonsyn
14) total_nonsyn
15) het-synonymous_variant
16) hom_alt-synonymous_variant
17) total_syn

# Usage: ./count_dSNP_per_individual.py [vcf_fp] [vep_fp] [bad_mut_del_vs_tol_fp] [out_dir] [out_prefix]
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
# VeP report generated for the same file
vep_fp = sys.argv[2]
# BAD_Mutations report with "Deleterious" vs "Tolerated" annotations
bad_mut_del_vs_tol_fp = sys.argv[3]
# Output directory
out_dir = os.path.expanduser(sys.argv[4])
# Output file prefix
out_prefix = sys.argv[5]

# Note: not very pretty but works and gets the info needed!
# Read VeP report
vep_list = []
with open(vep_fp, 'r') as f:
    for line in f:
        if line.startswith('##'):
            continue
        elif line.startswith('#Uploaded'):
            tmp = line.strip().split('\t')
            newline = ['id'] + tmp
            vep_list.append(newline)
        else:
            tmp = line.strip().split('\t')
            tmp_vid = tmp[0].split('_')
            var_id = '_'.join([tmp_vid[0], tmp_vid[1], tmp_vid[2]])
            newline = [var_id] + tmp
            vep_list.append(newline)
            
vep_df = pd.DataFrame(vep_list)
vep_df = vep_df.T.set_index(0).T
# Select only needed columns
subset_vep_df = vep_df[['id', 'Consequence']]

# Read bad_mut del vs tol
del_vs_tol_df = pd.read_table(bad_mut_del_vs_tol_fp, delimiter='\t', header=0)
# Add a unique ID that can be used with other datasets
uniqid_del_vs_tol = []
for vid in del_vs_tol_df['VariantID']:
    tmp = vid.split('_')
    newid = '_'.join([tmp[0], tmp[1], tmp[2]])
    uniqid_del_vs_tol.append(newid)

uniqid_del_vs_tol_df = pd.DataFrame(uniqid_del_vs_tol, columns=['id'])
new_del_vs_tol_df = pd.concat([uniqid_del_vs_tol_df, del_vs_tol_df], axis=1)
# Select only needed columns
subset_del_vs_tol_df = new_del_vs_tol_df[['id', 'Del_vs_Tol']]

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
counts_summary = [['sample', 'num_snps', 'num_hom_alt', 'num_het', 'prop_hom_alt', 'prop_het']]
for column in vcf_df:
    if column != "id" and column != "chr" and column != "pos":
        counts = vcf_df[column].value_counts()
        num_snps = counts['hom_alt']+counts['het']
        counts_summary.append([column, num_snps, counts['hom_alt'], counts['het'], counts['hom_alt']/num_snps, counts['het']/num_snps])

# Convert to pandas dataframe and set first row as column names
counts_summary_df = pd.DataFrame(counts_summary)
counts_summary_df = counts_summary_df.T.set_index(0).T

# Convert sample genotypes to rows with sample names
vcf_samples = pd.DataFrame()

for column in vcf_df:
    if column != "id" and column != "chr" and column != "pos":
        # column is the sample name
        tmp_curr_sample = vcf_df[['id', 'chr', 'pos', column]]
        # For current sample only keep 'het' and 'hom_alt', exclude 'miss' and 'hom_ref'
        # And rename column
        private_var = tmp_curr_sample[(tmp_curr_sample[column] != 'miss') & (tmp_curr_sample[column] != 'hom_ref')].rename(columns={column: 'gt'})
        # Add column of name of current sample
        private_var.insert(0, 'sample', column)
        # Append to data frame
        vcf_samples = vcf_samples.append(private_var, ignore_index=True)

# Add "Deleterious" vs "Tolerated" info (from BAD_Mutations results)
# subset_del_vs_tol_df contains only "id" and "Del_vs_Tol"
# Example: chr1H_part2_155053393	Deleterious
sample_del_vs_tol = subset_del_vs_tol_df.merge(vcf_samples, how='left', on=['id'])
# Save table to file
sample_del_vs_tol.to_csv(out_dir + '/' + out_prefix + '_per_sample_table.del_vs_tol.txt', sep='\t', header=True, na_rep='NA', index=False)

# Count the following categories for each sample:
# 1) het del
# 2) het tol
# 3) hom_alt del
# 4) hom_alt tol
# Create dataframe that stores count of each non-zero class in the new column 'counts'
counts_del_vs_tol = sample_del_vs_tol.groupby(['sample', 'gt', 'Del_vs_Tol']).size().reset_index(name='counts')
# Use pivot_table to get desired dataframe with counts for both existing and non-existing classes
tmp_final_counts_del_vs_tol = pd.pivot_table(counts_del_vs_tol, index=['sample', 'gt', 'Del_vs_Tol'], values='counts', fill_value=0, dropna=False, aggfunc=np.sum).reset_index()
# Create final counts table column names
tmp_final_counts_del_vs_tol['gt_del_vs_tol'] = tmp_final_counts_del_vs_tol['gt'] + '-' + tmp_final_counts_del_vs_tol['Del_vs_Tol']
final_counts_del_vs_tol_df = pd.pivot(tmp_final_counts_del_vs_tol, index=['sample'], columns=['gt_del_vs_tol'], values='counts').reset_index()
# Add columns for per sample total del and total tolerated based on summing het and hom_alt
final_counts_del_vs_tol_df['total_del'] = final_counts_del_vs_tol_df['het-Deleterious'] + final_counts_del_vs_tol_df['hom_alt-Deleterious']
final_counts_del_vs_tol_df['total_tol'] = final_counts_del_vs_tol_df['het-Tolerated'] + final_counts_del_vs_tol_df['hom_alt-Tolerated']

# Pull out synonymous variants from vep report
vep_syn = subset_vep_df[subset_vep_df['Consequence'] == 'synonymous_variant']
# Add synonymous label
sample_syn = vcf_samples.merge(vep_syn, how='inner', on=['id']).drop_duplicates()

# Count the following categories for each sample:
# 1) het syn
# 2) hom_alt syn
# Create dataframe that stores count of each non-zero class in the new column 'counts'
counts_syn = sample_syn.groupby(['sample', 'gt', 'Consequence']).size().reset_index(name='syn_counts')
# Use pivot_table to get desired dataframe with counts for both existing and non-existing classes
tmp_final_counts_syn = pd.pivot_table(counts_syn, index=['sample', 'gt', 'Consequence'], values='syn_counts', fill_value=0, dropna=False, aggfunc=np.sum).reset_index()
# Create final counts table column names
tmp_final_counts_syn['gt_vep_conseq'] = tmp_final_counts_syn['gt'] + '-' + tmp_final_counts_syn['Consequence']
final_counts_syn_df = pd.pivot(tmp_final_counts_syn, index=['sample'], columns=['gt_vep_conseq'], values='syn_counts').reset_index()
# Add columns for per sample total syn
final_counts_syn_df['total_syn'] = final_counts_syn_df['het-synonymous_variant'] + final_counts_syn_df['hom_alt-synonymous_variant']

# Pull out nonsynonymous variants which include consequences:
# missense_variant, start_lost, stop_gained, and stop_lost
vep_nonsyn = subset_vep_df[(subset_vep_df["Consequence"] == "missense_variant") | (subset_vep_df["Consequence"] == "start_lost") | (subset_vep_df["Consequence"] == "stop_gained") | (subset_vep_df["Consequence"] == "stop_lost")]
# Add column label for "nonsynonymous"
vep_nonsyn.insert(2, 'var_cat', 'nonsyn')
vep_nonsyn
# Add nonsynonymous label
sample_nonsyn = vcf_samples.merge(vep_nonsyn, how='inner', on=['id']).drop_duplicates()

# Count the following categories for each sample:
# 1) het nonsyn
# 2) hom_alt nonsyn
# Create dataframe that stores count of each non-zero class in the new column 'counts'
counts_nonsyn = sample_nonsyn.groupby(['sample', 'gt', 'var_cat']).size().reset_index(name='nonsyn_counts')
# Use pivot_table to get desired dataframe with counts for both existing and non-existing classes
tmp_final_counts_nonsyn = pd.pivot_table(counts_nonsyn, index=['sample', 'gt', 'var_cat'], values='nonsyn_counts', fill_value=0, dropna=False, aggfunc=np.sum).reset_index()
# Create final counts table column names
tmp_final_counts_nonsyn['gt_vep_conseq'] = tmp_final_counts_nonsyn['gt'] + '-' + tmp_final_counts_nonsyn['var_cat']
final_counts_nonsyn_df = pd.pivot(tmp_final_counts_nonsyn, index=['sample'], columns=['gt_vep_conseq'], values='nonsyn_counts').reset_index()
# Add columns for per sample total syn
final_counts_nonsyn_df['total_nonsyn'] = final_counts_nonsyn_df['het-nonsyn'] + final_counts_nonsyn_df['hom_alt-nonsyn']

# Combine all counts tables by sample
counts_vep_conseq = final_counts_nonsyn_df.merge(final_counts_syn_df, how='inner', on='sample')
counts_del_vep = final_counts_del_vs_tol_df.merge(counts_vep_conseq, how='inner', on='sample')
final_counts_combined = counts_summary_df.merge(counts_del_vep, how='inner', on='sample')

# Save final counts to table
final_counts_combined.to_csv(out_dir + '/' + out_prefix + '_per_sample_counts.del_vs_tol.nonsyn_vs_syn.txt', sep='\t', header=True, na_rep='NA', index=False)
