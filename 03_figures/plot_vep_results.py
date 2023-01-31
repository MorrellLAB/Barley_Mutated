#!/usr/bin/env python3
"""Code used to generate plot of vep results for mutated, hybrid rare, and hybrid common
1bp SNPs and INDEls."""

from cyvcf2 import VCF
import pandas as pd
from pandas.api.types import CategoricalDtype
import matplotlib.pyplot as plt
import numpy as np
import plotly.io as pio
import plotly.express as px

# Filepaths to VeP results (duplicate lines removed)
mut_snps_fp = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/HC_LC_gff_SNPs_private_all_samples/mut8_and_3mut10xGenomics.SNPs.private.txt"
mut_indels_fp = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/HC_LC_gff_INDELs_private_all_samples/mut8_and_3mut10xGenomics.INDELs.private.txt"

hyb_rare_snps_fp = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/HC_LC_gff_SNPs_hybrid13_rare/hybrid13.SNPs.rare.txt"
hyb_rare_indels_fp = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/HC_LC_gff_INDELs_hybrid13_rare/hybrid13.INDELs.rare.txt"

hyb_comm_snps_fp = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/HC_LC_gff_SNPs_hybrid13_common/hybrid13.SNPs.common.txt"
hyb_comm_indels_fp = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/HC_LC_gff_INDELs_hybrid13_common/hybrid13.INDELs.common.txt"

#--------
# Prepare unique consequences for each vep file
# Easier to do this on bash side:
# grep -v "#" mut8_and_3mut10xGenomics.SNPs.private.txt | cut -f 7 | sort -uV | tr '\n' ', ' | sed 's/\(.*\),/\1\n/' | sed 's/,/, /g' | sed 's/[^[:space:],]\+/"&"/g'
# Mutated
mut_snp_conseq_cat = ["3_prime_UTR_variant", "5_prime_UTR_variant", "downstream_gene_variant",
                        "intergenic_variant", "intron_variant", "missense_variant",
                        "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant",
                        "start_lost", "stop_gained", "stop_lost", "synonymous_variant", "upstream_gene_variant"]

mut_indel_conseq_cat = ["3_prime_UTR_variant", "5_prime_UTR_variant", "downstream_gene_variant",
                        "frameshift_variant", "intergenic_variant", "intron_variant", "splice_region_variant",
                        "upstream_gene_variant"]

# Hybrid rare
rare_snp_conseq_cat = ["3_prime_UTR_variant", "5_prime_UTR_variant", "coding_sequence_variant",
                       "downstream_gene_variant", "incomplete_terminal_codon_variant", "intergenic_variant",
                       "intron_variant", "missense_variant", "splice_acceptor_variant", "splice_donor_variant",
                       "splice_region_variant", "start_lost", "stop_gained", "stop_lost", "stop_retained_variant",
                       "synonymous_variant", "upstream_gene_variant"]

rare_indel_conseq_cat = ["3_prime_UTR_variant", "5_prime_UTR_variant", "downstream_gene_variant",
                         "frameshift_variant", "intergenic_variant", "intron_variant", "splice_acceptor_variant",
                         "splice_donor_variant", "splice_region_variant", "start_lost", "start_lost",
                         "start_retained_variant", "stop_gained", "stop_lost", "stop_retained_variant",
                         "upstream_gene_variant"]

# Hybrid common
comm_snp_conseq_cat = ["3_prime_UTR_variant", "5_prime_UTR_variant", "coding_sequence_variant",
                       "downstream_gene_variant", "incomplete_terminal_codon_variant", "intergenic_variant",
                       "intron_variant", "missense_variant", "splice_acceptor_variant", "splice_donor_variant",
                       "splice_region_variant", "start_lost", "stop_gained", "stop_lost", "stop_retained_variant",
                       "synonymous_variant", "upstream_gene_variant"]

comm_indel_conseq_cat = ["3_prime_UTR_variant", "5_prime_UTR_variant", "downstream_gene_variant",
                         "frameshift_variant", "intergenic_variant", "intron_variant", "splice_acceptor_variant",
                         "splice_donor_variant", "splice_region_variant", "start_lost", "start_lost",
                         "start_retained_variant", "stop_gained", "stop_lost", "stop_retained_variant",
                         "upstream_gene_variant"]

# All VeP categories
def prep_vep_report(vep_fp, var_consequences, data_name, var_type):
    """Prepare data for plotting"""
    # Prepare dict of consequence types
    conseq_dict = {}
    for i in var_consequences:
        conseq_dict[i] = 0
    # Prepare dict of consequence impact
    impact_dict = {}
    # Keep track of total (for calculating proportions later)
    total = 0
    
    with open(vep_fp, 'r') as f:
        for i in f:
            if not i.startswith("#"):
                v = i.strip().split()
                impact = v[13].split(';')[0].split('=')[1].title()
                c = v[6].split(',')
                # Increment by 1 for consequence type
                for cidx in range(0, len(c)):
                    total += 1
                    conseq_dict[c[cidx]] += 1
                    if c[cidx] not in impact_dict.keys():
                        impact_dict[c[cidx]] = impact
    # Calculate proportions
    var_prop = []
    for k in conseq_dict.keys():
        prop = conseq_dict[k]/total
        percent = round(prop * 100, 2)
        if prop < 0.001:
            # Add dataset label to bar label if bar will be too small to show up in plot
            var_prop.append([data_name, var_type, k, impact_dict[k], conseq_dict[k], total,
                     prop, percent,''.join([str(percent), '% (', str(conseq_dict[k]), ') - ', data_name])])
        else:
            var_prop.append([data_name, var_type, k, impact_dict[k], conseq_dict[k], total,
                             prop, percent,''.join([str(percent), '% (', str(conseq_dict[k]), ')'])])
    
    return var_prop

# Load and prepare vep reports
mut_snp_prop = prep_vep_report(mut_snps_fp, mut_snp_conseq_cat, 'Mutated', 'snp')
mut_indel_prop = prep_vep_report(mut_indels_fp, mut_indel_conseq_cat, 'Mutated', 'indel')
rare_snp_prop = prep_vep_report(hyb_rare_snps_fp, rare_snp_conseq_cat, 'Rare', 'snp')
rare_indel_prop = prep_vep_report(hyb_rare_indels_fp, rare_indel_conseq_cat, 'Rare', 'indel')
comm_snp_prop = prep_vep_report(hyb_comm_snps_fp, comm_snp_conseq_cat, 'Common', 'snp')
comm_indel_prop = prep_vep_report(hyb_comm_indels_fp, comm_indel_conseq_cat, 'Common', 'indel')

# Combine list of lists
combined = mut_snp_prop + mut_indel_prop + rare_snp_prop + rare_indel_prop + comm_snp_prop + comm_indel_prop

# Create custom sort for dataset
dataset_order = CategoricalDtype(['Mutated', 'Rare', 'Common'], ordered=True)
impact_order = CategoricalDtype(['Modifier', 'Low', 'Moderate', 'High'], ordered=True)

# Convert to pandas dataframe
combined_df = pd.DataFrame(combined, columns=['Dataset', 'Var_Type', 'Consequence', 'Impact', 'Count',
                                              'Total', 'Proportion', 'Percent', 'Bar_Label'])
combined_df['Dataset'] = combined_df['Dataset'].astype(dataset_order)
# Custom sort by impact categories
combined_df['Impact'] = combined_df['Impact'].astype(impact_order)
combined_sorted_df = combined_df.sort_values(['Impact', 'Dataset'])

# Plot all VeP categories - SNPs
snps_df = combined_sorted_df[combined_sorted_df['Var_Type'] == 'snp']

pio.templates.default = "plotly_white"

fig = px.bar(snps_df,
             x="Percent", y="Consequence", color="Dataset", text="Bar_Label",
             color_discrete_sequence=["#d22525", "#cee3f8", "#6898d6"])

fig.update_traces(textfont_size=10, textposition="outside")
fig.update_layout(barmode='group',
                  height = 1000,
                  uniformtext_minsize=10, uniformtext_mode='show',
                  xaxis_title='% (Count)',
                  yaxis_title='Consequence',
                  legend_title="",
                  legend=dict(
                      yanchor="top", y=0.99,
                      xanchor="right", x=0.99),
                 xaxis=dict(range=[0, max(snps_df['Percent'])+max(snps_df['Percent']*0.16)]))
fig.update_xaxes(automargin=True)
fig.show()

# Plot all VeP categories - INDELs
indels_df = combined_sorted_df[combined_sorted_df['Var_Type'] == 'indel']

pio.templates.default = "plotly_white"

fig = px.bar(indels_df,
             x="Percent", y="Consequence", color="Dataset", text="Bar_Label",
             color_discrete_sequence=["#d22525", "#cee3f8", "#6898d6"])

fig.update_traces(textfont_size=10, textposition="outside")
fig.update_layout(barmode='group',
                  height = 1000,
                  uniformtext_minsize=10, uniformtext_mode='show',
                  xaxis_title='% (Count)',
                  yaxis_title='Consequence',
                  legend_title="",
                  legend=dict(
                      yanchor="top", y=0.99,
                      xanchor="right", x=0.99),
                 xaxis=dict(range=[0, max(snps_df['Percent'])+max(snps_df['Percent']*0.16)]))
fig.update_xaxes(automargin=True)
fig.show()

# Plot VeP categories with some binned
# Bin the following based on impact categories
#   splice sites, start and stop, and 5' and 3'
# Work with list of lists already created above
def bin_vep_cat(dat_list):
    tmp_bins = {
    '5_prime_UTR_variant,3_prime_UTR_variant': [],
    'splice_acceptor_variant,splice_donor_variant': [],
    'start_lost,stop_gained,stop_lost,start_retained_variant,stop_retained_variant': []
    }
    
    output_bins = []
    
    for i in dat_list:
        # For the following categories, bin together
        if i[2] == '3_prime_UTR_variant' or i[2] == '5_prime_UTR_variant':
            tmp_bins['5_prime_UTR_variant,3_prime_UTR_variant'].append(i)
        elif i[2] == 'splice_acceptor_variant' or i[2] == 'splice_donor_variant':
            tmp_bins['splice_acceptor_variant,splice_donor_variant'].append(i)
        elif i[2] == 'start_lost' or i[2] == 'stop_gained' or i[2] == 'stop_lost' or i[2] == 'stop_retained_variant' or i[2] == 'start_retained_variant':
            tmp_bins['start_lost,stop_gained,stop_lost,start_retained_variant,stop_retained_variant'].append(i)
        else:
            output_bins.append(i)
    
    new_bins = {}
    
    for key in tmp_bins.keys():
        tmp_sum = 0
        if len(tmp_bins[key]) != 0:
            for i in tmp_bins[key]:
                tmp_sum += i[4]
            # Combine info
            new_prop = tmp_sum/tmp_bins[key][0][5]
            new_percent = round(new_prop * 100, 2)
            if new_prop < 0.001:
                # Add dataset label to bar label
                new_bins[key] = [tmp_bins[key][0][0], tmp_bins[key][0][1], key, tmp_bins[key][0][3],
                                tmp_sum, tmp_bins[key][0][5], new_prop, new_percent,
                                ''.join([str(new_percent), '% (', str(tmp_sum), ') - ', tmp_bins[key][0][0]])]
            else:
                new_bins[key] = [tmp_bins[key][0][0], tmp_bins[key][0][1], key, tmp_bins[key][0][3],
                                tmp_sum, tmp_bins[key][0][5], new_prop, new_percent,
                                ''.join([str(new_percent), '% (', str(tmp_sum), ')'])]

    for key in new_bins.keys():
        output_bins.append(new_bins[key])
    return output_bins

bin_mut_snp_prop = bin_vep_cat(mut_snp_prop)
bin_mut_indel_prop = bin_vep_cat(mut_indel_prop)

bin_rare_snp_prop = bin_vep_cat(rare_snp_prop)
bin_rare_indel_prop = bin_vep_cat(rare_indel_prop)

bin_comm_snp_prop = bin_vep_cat(comm_snp_prop)
bin_comm_indel_prop = bin_vep_cat(comm_indel_prop)

# For categories that only appear in 2 out of 3 datasets
# add in placeholder row so plot appears in the correct order
bin_mut_snp_prop.append(['Mutated', 'snp', 'incomplete_terminal_codon_variant', 'Low', 0, 43706, 0, 0, '0% (0) - Mutated'])
bin_mut_snp_prop.append(['Mutated', 'snp', 'coding_sequence_variant', 'Low', 0, 43706, 0, 0, '0% (0) - Mutated'])

bin_mut_indel_prop.append(['Mutated', 'indel', 'start_lost,stop_gained,stop_lost,start_retained_variant,stop_retained_variant', 'High', 0, 11832, 0, 0, '0% (0) - Mutated'])
bin_mut_indel_prop.append(['Mutated', 'indel', 'splice_acceptor_variant,splice_donor_variant', 'High', 0, 11832, 0, 0, '0% (0) - Mutated'])

bin_combined_prop = bin_mut_snp_prop + bin_mut_indel_prop + bin_rare_snp_prop + bin_rare_indel_prop + bin_comm_snp_prop + bin_comm_indel_prop

# Create custom sort based on impact and consequence
conseq_order = CategoricalDtype([
    'intergenic_variant', 'downstream_gene_variant', 'upstream_gene_variant', 'intron_variant', '5_prime_UTR_variant,3_prime_UTR_variant',
    'synonymous_variant', 'splice_region_variant', 'coding_sequence_variant', 'incomplete_terminal_codon_variant',
    'missense_variant',
    'frameshift_variant', 'splice_acceptor_variant,splice_donor_variant', 'start_lost,stop_gained,stop_lost,start_retained_variant,stop_retained_variant'],
    ordered=True)


df_mapping = pd.DataFrame({
    'conseq_order': [
    'intergenic_variant', 'downstream_gene_variant', 'upstream_gene_variant', 'intron_variant', '5_prime_UTR_variant,3_prime_UTR_variant',
    'synonymous_variant', 'splice_region_variant', 'coding_sequence_variant', 'incomplete_terminal_codon_variant',
    'missense_variant',
    'frameshift_variant', 'splice_acceptor_variant,splice_donor_variant', 'start_lost,stop_gained,stop_lost,start_retained_variant,stop_retained_variant'],
})
df_sort_conseq_mapping = df_mapping.reset_index().set_index('conseq_order')

bin_combined_df = pd.DataFrame(bin_combined_prop, columns=['Dataset', 'Var_Type', 'Consequence', 'Impact', 'Count',
                                              'Total', 'Proportion', 'Percent', 'Bar_Label'])
bin_combined_df['Dataset'] = bin_combined_df['Dataset'].astype(dataset_order)
# Custom sort by impact categories
bin_combined_df['conseq_order'] = bin_combined_df['Consequence'].map(df_sort_conseq_mapping['index'])
bin_combined_sorted_df = bin_combined_df.sort_values(['conseq_order', 'Dataset'])

# SNPs
bin_snps_df = bin_combined_sorted_df[bin_combined_sorted_df['Var_Type'] == 'snp']

pio.templates.default = "plotly_white"

fig = px.bar(bin_snps_df,
             x="Percent", y="Consequence", color="Dataset", text="Bar_Label",
             color_discrete_sequence=["#d22525", "#cee3f8", "#6898d6"])

fig.update_traces(textfont_size=12, textposition="outside")
fig.update_layout(barmode='group',
                  width = 1400,
                  height = 1000,
                  uniformtext_minsize=10, uniformtext_mode='show',
                  xaxis_title='% (Count)',
                  yaxis_title='',
                  legend_title="",
                  legend=dict(
                      yanchor="top", y=0.99,
                      xanchor="right", x=0.99),
                 #xaxis=dict(range=[0, max(bin_snps_df['Percent'])+max(bin_snps_df['Percent']*0.16)])
                 xaxis=dict(range=[0, 50]))
fig.update_xaxes(automargin=True)
fig.show()
# Save image to file
fig.write_image("/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures/vep_summary-snps.png")
fig.write_image("/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures/vep_summary-snps.svg")


# Indels
bin_indels_df = bin_combined_sorted_df[bin_combined_sorted_df['Var_Type'] == 'indel']

pio.templates.default = "plotly_white"

fig = px.bar(bin_indels_df,
             x="Percent", y="Consequence", color="Dataset", text="Bar_Label",
             color_discrete_sequence=["#d22525", "#cee3f8", "#6898d6"])

fig.update_traces(textfont_size=12, textposition="outside")
fig.update_layout(barmode='group',
                  width = 1400,
                  height = 700,
                  uniformtext_minsize=10, uniformtext_mode='show',
                  xaxis_title='% (Count)',
                  yaxis_title='',
                  legend_title="",
                  legend=dict(
                      yanchor="top", y=0.99,
                      xanchor="right", x=0.99),
                 #xaxis=dict(range=[0, max(bin_snps_df['Percent'])+max(bin_snps_df['Percent']*0.16)])
                 xaxis=dict(range=[0, 50]))
fig.update_xaxes(automargin=True)
fig.show()
# Save image to file
fig.write_image("/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures/vep_summary-indels.png")
fig.write_image("/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures/vep_summary-indels.svg")
