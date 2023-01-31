#!/usr/bin/env python3
"""Code used to generate sunburst plot of mut, hybrid rare, and hybrid common
1bp SNPs and INDEls."""

from cyvcf2 import VCF
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import plotly.io as pio
import plotly.express as px

# List of VCF files
# Mut lines private variants (1bp)
mut_snps_vcf_list_fp = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/split_by_sample_SNPs_private/mut_snps_private_list.txt"
mut_indels_vcf_list_fp = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/split_by_sample_INDELs_private/mut_indels_private_list.txt"

# Hybrid rare variants (1bp)
hyb_rare_snps_vcf_list_fp = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/split_by_sample_SNPs_rare/hybrid13_snps_rare_list.txt"
hyb_rare_indels_vcf_list_fp = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/split_by_sample_INDELs_rare/hybrid13_indels_rare_list.txt"

# Hybrid common variants (1bp)
hyb_common_snps_vcf_list_fp = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/split_by_sample_SNPs_common/hybrid13_snps_common_list.txt"
hyb_common_indels_vcf_list_fp = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/split_by_sample_INDELs_common/hybrid13_indels_common_list.txt"

# SNPs and INDELs that overlap GFF
mut_snps_gff_fp = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.SNPs.private.gffOverlap.vcf"
mut_indels_gff_fp = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.INDELs.private.gffOverlap.vcf"
hyb_snps_rare_gff_fp = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/hybrid13.SNPs.rare.gffOverlap.vcf"
hyb_indels_rare_gff_fp = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/hybrid13.INDELs.rare.gffOverlap.vcf"
hyb_snps_comm_gff_fp = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/hybrid13.SNPs.common.gffOverlap.vcf"
hyb_indels_comm_gff_fp = "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/hybrid13.INDELs.common.gffOverlap.vcf"

def prep_gff_overlap_info(vcf_fp, curr_var_type):
    """Prepare dictionary to lookup variants that overlap GFF."""
    curr_vcf = VCF(vcf_fp)
    vcf_dict = {}
    for v in curr_vcf:
        vcf_dict['_'.join([v.CHROM, str(v.POS)])] = curr_var_type
    return vcf_dict


def prep_sample_vcf_info(vcf_list_fp, gff_overlap_dict, curr_var_type):
    """Prepare list of lists containing info about sample and hom/het.
    Each VCF file in list should contain a single sample and no missing data."""
    var_list = []
    # Read in list of VCF file paths
    vcf_list = open(vcf_list_fp, 'r').readlines()
    # Prepare info
    for vfile in vcf_list:
        curr_vcf = VCF(vfile.strip())
        curr_sample = curr_vcf.samples[0]
        nhet = 0
        nhom = 0
        pot_dam_het = 0
        pot_dam_hom = 0
        not_dam_het = 0
        not_dam_hom = 0
        tmp_var_list = []
        for variant in curr_vcf:
            if variant.genotypes[0][0] == variant.genotypes[0][1]:
                h = 'hom'
            elif variant.genotypes[0][0] != variant.genotypes[0][1]:
                h = 'het'
            #var_list.append([curr_sample, curr_var_type, variant.CHROM, variant.POS, variant.genotypes[0], h])
            tmp_var_list.append([curr_sample, curr_var_type, variant.CHROM, variant.POS, variant.genotypes[0], h])
        # Count number of hom and het for each sample
        for g in tmp_var_list:
            tmp_chr_pos = '_'.join([g[2], str(g[3])])
            if g[5] == 'het':
                nhet += 1
                if tmp_chr_pos in gff_overlap_dict.keys():
                    pot_dam_het += 1
                elif tmp_chr_pos not in gff_overlap_dict.keys():
                    not_dam_het += 1
            elif g[5] == 'hom':
                nhom += 1
                if tmp_chr_pos in gff_overlap_dict.keys():
                    pot_dam_hom += 1
                elif tmp_chr_pos not in gff_overlap_dict.keys():
                    not_dam_hom += 1
        #var_list.append([curr_sample, curr_var_type, 'total', 'het', nhet, pot_dam_het, not_dam_het])
        #var_list.append([curr_sample, curr_var_type, 'total', 'hom', nhom, pot_dam_hom, not_dam_hom])
        var_list.append([curr_sample, curr_var_type, 'total', 'het', nhet, 'Potentially Damaging', pot_dam_het])
        var_list.append([curr_sample, curr_var_type, 'total', 'het', nhet, 'Not Damaging', not_dam_het])
        var_list.append([curr_sample, curr_var_type, 'total', 'hom', nhom, 'Potentially Damaging', pot_dam_hom])
        var_list.append([curr_sample, curr_var_type, 'total', 'hom', nhom, 'Not Damaging', not_dam_hom])
    return var_list

# Prepare dict to store variants that overlap with GFF
mut_snps_gff = prep_gff_overlap_info(mut_snps_gff_fp, 'snp')
mut_indels_gff = prep_gff_overlap_info(mut_indels_gff_fp, 'indel')

hyb_snps_rare_gff = prep_gff_overlap_info(hyb_snps_rare_gff_fp, 'snp')
hyb_indels_rare_gff = prep_gff_overlap_info(hyb_indels_rare_gff_fp, 'indel')

hyb_snps_comm_gff = prep_gff_overlap_info(hyb_snps_comm_gff_fp, 'snp')
hyb_indels_comm_gff = prep_gff_overlap_info(hyb_indels_comm_gff_fp, 'indel')

# Prepare dict to store variants that overlap with GFF
mut_snps_list = prep_sample_vcf_info(mut_snps_vcf_list_fp, mut_snps_gff, 'snp')
mut_indels_list = prep_sample_vcf_info(mut_indels_vcf_list_fp, mut_indels_gff, 'indel')

hyb_rare_snps_list = prep_sample_vcf_info(hyb_rare_snps_vcf_list_fp, hyb_snps_rare_gff, 'snp')
hyb_rare_indels_list = prep_sample_vcf_info(hyb_rare_indels_vcf_list_fp, hyb_indels_rare_gff, 'indel')

hyb_common_snps_list = prep_sample_vcf_info(hyb_common_snps_vcf_list_fp, hyb_snps_comm_gff, 'snp')
hyb_common_indels_list = prep_sample_vcf_info(hyb_common_indels_vcf_list_fp, hyb_indels_comm_gff, 'indel')

# Combine lists
mut = mut_snps_list + mut_indels_list
hyb_rare = hyb_rare_snps_list + hyb_rare_indels_list
hyb_common = hyb_common_snps_list + hyb_common_indels_list

# Convert to pandas DataFrame for plotting
# Create Pandas dataframe for plotting
mut_df = pd.DataFrame(mut, columns=['Sample', 'Var_Type', 'Total', 'GT_Type', 'Count', 'Pot_Dam_Cat', 'Num_Pot_Dam'])
hyb_rare_df = pd.DataFrame(hyb_rare, columns=['Sample', 'Var_Type', 'Total', 'GT_Type', 'Count', 'Pot_Dam_Cat', 'Num_Pot_Dam'])
hyb_common_df = pd.DataFrame(hyb_common, columns=['Sample', 'Var_Type', 'Total', 'GT_Type', 'Count', 'Pot_Dam_Cat', 'Num_Pot_Dam'])

# Sunburst plot
def plot_tvs(df, out_fp):
    fig = px.sunburst(df, path=['Total', 'Var_Type', 'Sample'], values='Count', color='Var_Type',
                 color_discrete_sequence=['white', '#d3bec2', '#e4eaef'])
    fig.update_traces(textinfo="label+value+percent parent")
    fig.show()
    fig.write_image(out_fp)

# Generate sunburst plot and save to file
# Saving to SVG allows easy modifying of font size and positioning in OmniGraffle
plot_tvs(mut_df, "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures/sunburst-mut.png")
plot_tvs(mut_df, "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures/sunburst-mut.svg")

plot_tvs(hyb_rare_df, "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures/sunburst-hyb_rare.png")
plot_tvs(hyb_rare_df, "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures/sunburst-hyb_rare.svg")

plot_tvs(hyb_common_df, "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures/sunburst-hyb_common.png")
plot_tvs(hyb_common_df, "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures/sunburst-hyb_common.svg")

# Potentially damaging SNPs and 1 bp INDELs
def prep_dam_info(df):
    """Reformat dataframe and calculate proportions for plotting."""
    # Get number of potentially damaging vs not for each sample's SNPs and INDELs
    dam_df = df.groupby(['Sample', 'Var_Type', 'Pot_Dam_Cat'])['Num_Pot_Dam'].sum().reset_index()
    # Calculate proportion
    tmp_prop = df.groupby(['Sample', 'Var_Type', 'Pot_Dam_Cat'])['Num_Pot_Dam'].sum() / df.groupby(['Sample', 'Var_Type'])['Num_Pot_Dam'].sum()
    # Add proportion column
    dam_df['Proportion'] = tmp_prop.reset_index()['Num_Pot_Dam']
    # Add percentage (for labeling bars)
    dam_df['Percent'] = round(dam_df['Proportion'] * 100, 1)
    # Add new label
    dam_df['New_Label'] = dam_df['Num_Pot_Dam'].astype('str') + ' (' + dam_df['Percent'].astype('str') + '%)'
    return dam_df


def plot_pot_dam_snps(df, out_fp):
    snps_df = df[df['Var_Type'] == 'snp'].sort_values('Num_Pot_Dam')
    
    pio.templates.default = "plotly_white"
    
    fig = px.bar(snps_df,
                 x="Num_Pot_Dam", y="Sample", color="Pot_Dam_Cat", text="New_Label",
                 color_discrete_sequence=['#d3bec2', '#e4eaef'])

    fig.update_traces(textfont_size=12, textposition="outside")
    fig.update_layout(barmode='group',
                      uniformtext_minsize=10, uniformtext_mode='show',
                      xaxis_title='Number of SNPs',
                      yaxis_title='Sample',
                      legend_title="",
                      legend=dict(
                          yanchor="bottom", y=0.01,
                          xanchor="right", x=0.99),
                      xaxis=dict(range=[0, max(snps_df['Num_Pot_Dam'])+max(snps_df['Num_Pot_Dam']*0.16)]))
    fig.update_xaxes(automargin=True)
    fig.show()
    fig.write_image(out_fp)


def plot_pot_dam_indels(df, out_fp):
    indels_df = df[df['Var_Type'] == 'indel'].sort_values('Num_Pot_Dam')
    
    pio.templates.default = "plotly_white"
    
    fig = px.bar(indels_df,
                 x="Num_Pot_Dam", y="Sample", color="Pot_Dam_Cat", text="New_Label",
                 color_discrete_sequence=['#d3bec2', '#e4eaef'])

    fig.update_traces(textfont_size=12, textposition="outside")
    fig.update_layout(barmode='group',
                      uniformtext_minsize=10, uniformtext_mode='show',
                      xaxis_title='Number of INDELs',
                      yaxis_title='Sample',
                      legend_title="",
                      legend=dict(
                          yanchor="bottom", y=0.01,
                          xanchor="right", x=0.99),
                     xaxis=dict(range=[0, max(indels_df['Num_Pot_Dam'])+max(indels_df['Num_Pot_Dam']*0.16)]))
    fig.update_xaxes(automargin=True)
    fig.show()
    fig.write_image(out_fp)

# Prepare data
dam_mut_df = prep_dam_info(mut_df)
dam_hyb_rare_df = prep_dam_info(hyb_rare_df)
dam_hyb_comm_df = prep_dam_info(hyb_common_df)

# Generate grouped barplots per sample for potentially damaging SNPs and 1 bp INDELs
# Save to PNG and SVG
# Mutated
plot_pot_dam_snps(dam_mut_df,
                  "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures/potentially_damaging_snps-mut.png")
plot_pot_dam_snps(dam_mut_df,
                  "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures/potentially_damaging_snps-mut.svg")

plot_pot_dam_indels(dam_mut_df,
                   "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures/potentially_damaging_indels-mut.png")
plot_pot_dam_indels(dam_mut_df,
                   "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures/potentially_damaging_indels-mut.svg")

# Hybrid rare
plot_pot_dam_snps(dam_hyb_rare_df,
                 "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures/potentially_damaging_snps-hyb_rare.png")
plot_pot_dam_snps(dam_hyb_rare_df,
                 "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures/potentially_damaging_snps-hyb_rare.svg")

plot_pot_dam_indels(dam_hyb_rare_df,
                   "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures/potentially_damaging_indels-hyb_rare.png")
plot_pot_dam_indels(dam_hyb_rare_df,
                   "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures/potentially_damaging_indels-hyb_rare.svg")

# Hybrid common
plot_pot_dam_snps(dam_hyb_comm_df,
                 "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures/potentially_damaging_snps-hyb_comm.png")
plot_pot_dam_snps(dam_hyb_comm_df,
                 "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures/potentially_damaging_snps-hyb_comm.svg")

plot_pot_dam_indels(dam_hyb_comm_df,
                   "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures/potentially_damaging_indels-hyb_comm.png")
plot_pot_dam_indels(dam_hyb_comm_df,
                   "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/figures/potentially_damaging_indels-hyb_comm.svg")
