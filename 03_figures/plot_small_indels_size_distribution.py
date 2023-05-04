#!/usr/bin/env python3
"""Plot small indels size distribution

All mutated lines for "short read" data (Illumina WGS and 10x Genomics Longranger phased_variants.vcf file) compared to rare and common variants.
"""

from cyvcf2 import VCF
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Filepaths
mut_vcf_fp = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.INDELs.private.vcf.gz"
rare_vcf_fp = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/hybrid13.INDELs.rare.vcf.gz"
common_vcf_fp = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/hybrid13.INDELs.common.vcf.gz"
# Plot output directory
out_dir = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures"

#------------
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
    # Get the count in each size group
    count_df = indel_sizes_df.groupby('dataset')['indel_size'].value_counts().fillna(0).unstack(-1).transpose().reset_index()
    # Rename columns
    prop_df.rename(columns={dataset_name: "proportion"}, inplace=True)
    count_df.rename(columns={dataset_name: "count"}, inplace=True)
    # Merge df
    merged_df = count_df.merge(prop_df, how='left', on='indel_size')
    # Add dataset category
    merged_df['dataset'] = dataset_name
    return(merged_df)


mut_df = read_indels(mut_vcf_fp, "Sodium azide")
rare_df = read_indels(rare_vcf_fp, "Rare")
common_df = read_indels(common_vcf_fp, "Common")

# Combine dataframes
combined_df = pd.concat([mut_df, rare_df, common_df], ignore_index=True)
# Convert indel_size to numeric
combined_df['indel_size'] = combined_df['indel_size'].apply(pd.to_numeric)
# Add row for "0" as placeholder in plot
combined_df.loc[len(combined_df.index)] = [0, 0, 0, 'Sodium azide']
# Save to file
combined_df.to_csv(out_dir + '/small_indel_sizes_prop.txt', sep='\t', index=False)

# Print max indel sizes
print('Largest insertion - sodium azide: ', max(combined_df[combined_df['dataset'] == 'Sodium azide']['indel_size']))
print('Largest insertion - rare: ', max(combined_df[combined_df['dataset'] == 'Rare']['indel_size']))
print('Largest insertion - common: ', max(combined_df[combined_df['dataset'] == 'Common']['indel_size']))

print('Largest deletion - sodium azide: ', min(combined_df[combined_df['dataset'] == 'Sodium azide']['indel_size']))
print('Largest deletion - rare: ', min(combined_df[combined_df['dataset'] == 'Rare']['indel_size']))
print('Largest deletion - common: ', min(combined_df[combined_df['dataset'] == 'Common']['indel_size']))

# Select up to 20 bp sized indels
df_20bp_indels = combined_df[(combined_df['indel_size'] <= 20) & (combined_df['indel_size'] >= -20)]

# Output plot filepath
pfp = out_dir + "/small_indels_size_distribution.jpg"

# Generate figure
sns.set_theme(style="white", font_scale=2)
sns.set_style('ticks')

g = sns.catplot(data=df_20bp_indels,
               kind="bar",
               x="indel_size", y="proportion", hue="dataset",
               palette=["#d22525", "#cee3f8", "#6898d6"],
               height=8, width=0.95,
               aspect=2)

g.set_xticklabels(['-20', '', '', '', '',
                   '-15', '', '', '', '',
                   '-10', '', '', '', '',
                   '-5', '', '', '', '',
                   '0',
                   '', '', '', '',
                   '5', '', '', '', '',
                   '10', '', '', '', '',
                   '15', '', '', '', '',
                   '20'])
g.set(ylim=(0, max(df_20bp_indels['proportion']+(max(df_20bp_indels['proportion'])*0.11))))
g.set_axis_labels("Length (bp)", "Proportion of indels")
g.legend.set_title("")
sns.move_legend(g, "upper left", bbox_to_anchor=(0.1, 0.96), title=None, frameon=False)

# Save as high-resolution png
plt.savefig(pfp, dpi=300)
