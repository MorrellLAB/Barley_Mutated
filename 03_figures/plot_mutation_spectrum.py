#!/usr/bin/env python3
"""Code used to generate mutation spectrum plot"""

from cyvcf2 import VCF
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Filepaths
mut_vcf_fp = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.SNPs.private.vcf.gz"
hyb_rare_vcf_fp = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/hybrid13.SNPs.rare.vcf.gz"
hyb_common_vcf_fp = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/hybrid13.SNPs.common.vcf.gz"
# Plot output directory
out_dir = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures"

def count_mut(vcf_fp, data_grouping, clean_name):
    """Count transitions and transversions for plotting mutation spectrum."""
    # Transversions and reverse complement
    tv_rc = ('AC', 'TG', 'AT', 'TA', 'CA', 'GT', 'CG', 'GC')
    # Transitions and reverse complement
    ts_rc = ('AG', 'TC', 'CT', 'GA')

    total = 0
    # Ts and Tv bins includes reverse complement
    ag_ts = 0
    ct_ts = 0
    ac_tv = 0
    at_tv = 0
    ca_tv = 0
    cg_tv = 0
    # Count number in each bin
    for v in VCF(vcf_fp):
        total += 1
        m = ''.join([v.REF, v.ALT[0]])
        if m in ts_rc:
            # Transition
            if m == 'AG' or m == 'TC':
                ag_ts += 1
            elif m == 'CT' or m == 'GA':
                ct_ts += 1
        elif m in tv_rc:
            # Transversion
            if m == 'AC' or m == 'TG':
                ac_tv += 1
            elif m == 'AT' or m == 'TA':
                at_tv += 1
            elif m == 'CA' or m == 'GT':
                ca_tv += 1
            elif m == 'CG' or m == 'GC':
                cg_tv += 1
    # Prepare data for plotting
    mlist = []
    mlist.append([data_grouping, clean_name, 'Transitions', 'A->G*', ag_ts, ag_ts/total])
    mlist.append([data_grouping, clean_name, 'Transitions', 'C->T*', ct_ts, ct_ts/total])
    mlist.append([data_grouping, clean_name, 'Transversions', 'A->C*', ac_tv, ac_tv/total])
    mlist.append([data_grouping, clean_name, 'Transversions', 'A->T*', at_tv, at_tv/total])
    mlist.append([data_grouping, clean_name, 'Transversions', 'C->A*', ca_tv, ca_tv/total])
    mlist.append([data_grouping, clean_name, 'Transversions', 'C->G*', cg_tv, cg_tv/total])
    return mlist

# Count transitions and transversions for each dataset
mut_list = count_mut(mut_vcf_fp, '1mut', 'Sodium azide')
hyb_rare_list = count_mut(hyb_rare_vcf_fp, '2hyb_rare', 'Rare non-mutagenized')
hyb_common_list = count_mut(hyb_common_vcf_fp, '3hyb_common', 'Common non-mutagenized')

# Combine
ms_list = mut_list + hyb_rare_list + hyb_common_list

# Convert to pandas dataframe for plotting
mdf = pd.DataFrame(ms_list, columns=['Group', 'Clean_Name', 'Mut_Cat', 'Mut', 'Count', 'Proportion']).sort_values(by=['Group', 'Proportion'])
# Save to file
mdf.to_csv(out_dir + '/mutation_spectrum_counts_prop.txt', sep='\t', index=False)

# Output plot filepath
pfp = out_dir + "/mutation_spectrum.png"

# Set plot theme
sns.set_theme(style="white", font_scale=1.5)

# Draw a nested barplot by species and sex
g = sns.catplot(
    data=mdf, kind="bar",
    x="Mut", y="Proportion", hue="Clean_Name",
    palette=["#d22525", "#cee3f8", "#6898d6"], height=6,
    aspect=1.3
)
g.despine(left=True)
g.set_axis_labels("", "Proportion of SNPs")
g.legend.set_title("")
sns.move_legend(g, "upper left", bbox_to_anchor=(0.15, 0.96), title=None, frameon=False)

# Save as high-resolution png
plt.savefig(pfp, dpi=300)
