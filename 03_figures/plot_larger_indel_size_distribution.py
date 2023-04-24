#!/usr/bin/env python3
"""Plot size distribution of larger indels.

1) 10x Genomics larger deletions
2) Sniffles2 ONT mutated - indels
3) cuteSV ONT mutated - indels
"""

from cyvcf2 import VCF
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Filepaths
fp_mut_10xG_dels = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/de_novo_larger_svs/mut_3_lines_dels_merged.callable.noRefDiffs.private.supports.final.vcf"
fp_mut_10xG_larger_svs = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/de_novo_larger_svs/mut_3_lines_large_svs_merged.callable.noRefDiffs.vcf"

fp_mut_ont_sniff_indel = "/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/filtered/mut_ont.private.geSup3.callable.noRefDiffs.final.INDELs.vcf"

fp_mut_ont_cutesv_indel = "/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/cutesv_calls/filtered/mut_ont_cutesv.private.callable.noRefDiffs.final.INDELs.vcf"

#-------------------
def load_del_sizes(vcf_fp, dataset):
    dels_sizes = []
    for v in VCF(vcf_fp):
        dels_sizes.append([dataset, v.INFO.get('SVLEN')])
    return dels_sizes


def load_indel_sizes(vcf_fp, dataset):
    indel_sizes = []
    for v in VCF(vcf_fp):
        indel_sizes.append([dataset, v.INFO.get('SVLEN')])
    return indel_sizes


mut_10xG_dels = load_del_sizes(fp_mut_10xG_dels, '10xGenomics')
mut_10xG_larger_svs = load_del_sizes(fp_mut_10xG_larger_svs, '10xGenomics')
mut_ont_sniff_indels = load_indel_sizes(fp_mut_ont_sniff_indel, 'ONT Sniffles2')
mut_ont_cutesv_indels = load_indel_sizes(fp_mut_ont_cutesv_indel, 'ONT cuteSV')

indels_list = mut_10xG_dels + mut_10xG_larger_svs + mut_ont_sniff_indels + mut_ont_cutesv_indels

indels_df = pd.DataFrame(indels_list, columns=['Dataset', 'indel_size'])

# Pull deletions only
dels_df = indels_df[indels_df["indel_size"] <= 0]
ins_df = indels_df[indels_df["indel_size"] >= 0]

# Quick glance at min and max sizes
print(max(dels_df['indel_size']))
print(min(dels_df['indel_size']))

print(max(dels_df[dels_df['Dataset'] == '10xGenomics']['indel_size']))
print(min(dels_df[dels_df['Dataset'] == '10xGenomics']['indel_size']))

print(max(dels_df[dels_df['Dataset'] == 'ONT Sniffles2']['indel_size']))
print(min(dels_df[dels_df['Dataset'] == 'ONT Sniffles2']['indel_size']))

print(max(dels_df[dels_df['Dataset'] == 'ONT cuteSV']['indel_size']))
print(min(dels_df[dels_df['Dataset'] == 'ONT cuteSV']['indel_size']))

print(max(ins_df['indel_size']))
print(min(ins_df['indel_size']))

print(max(ins_df[ins_df['Dataset'] == 'ONT Sniffles2']['indel_size']))
print(min(ins_df[ins_df['Dataset'] == 'ONT Sniffles2']['indel_size']))

print(max(ins_df[ins_df['Dataset'] == 'ONT cuteSV']['indel_size']))
print(min(ins_df[ins_df['Dataset'] == 'ONT cuteSV']['indel_size']))

# Figure of deletions size distribution
sns.set_theme(style="white", font_scale=1.2)

g = sns.displot(dels_df, x="indel_size", hue="Dataset",
                hue_order=["ONT cuteSV", "ONT Sniffles2", "10xGenomics"],
                kind="kde", fill=True, aspect=1.5, clip=(0, min(dels_df['indel_size'])))

g.set(xlabel='Size of deletions (bp)')
sns.move_legend(g, "upper left", bbox_to_anchor=(0.2, 0.9))

plt.savefig("/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures/mut3_larger_dels_size_distribution.jpg",
            format="jpg", dpi=300)
plt.show()

# Figure of insertions size distribution
sns.set_theme(style="white", font_scale=1.2)

i = sns.displot(ins_df, x="indel_size", hue="Dataset",
                hue_order=["ONT cuteSV", "ONT Sniffles2"],
                kind="kde", fill=True, aspect=1.5, clip=(0, max(ins_df['indel_size'])))

i.set(xlabel='Size of insertions (bp)')
sns.move_legend(i, "upper right", bbox_to_anchor=(0.8, 0.9))
plt.savefig("/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures/mut3_larger_ins_size_distribution.jpg",
            format="jpg", dpi=300)
plt.show()

