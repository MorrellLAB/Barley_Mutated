#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=12gb
#SBATCH --tmp=6gb
#SBATCH -t 03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t,amd512,amd2tb
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Intersect SNP and indels called by 10x Genonmics Longranger to see if
#   there are sites where the call differs

# Dependencies
module load bedtools/2.29.2

# User provided input arguments
VCF_INVARIANT="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/Genotype_GVCFs/filtered/morex-sample2.invariant_sites.filtDP-RGQ-QUAL.vcf.gz"
# 10x Genomics filterd phased variants SNPs and INDELs
VCF_10X_SNPs="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/diffs_from_ref/morex-sample2_phased_variants-snps.DPfilt.noRepeatOverlap.noRefNs.vcf.gz"
VCF_10X_INDELs="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/diffs_from_ref/morex-sample2_phased_variants-indels.DPfilt.noRepeatOverlap.noRefNs.vcf.gz"
OUT_DIR="/scratch.global/liux1299/temp_morex-sample2"
OUT_PREFIX="morex-sample2"

#-----------------
mkdir -p ${OUT_DIR}

# Intersect vcf files
bedtools intersect -sorted -wa -a ${VCF_10X_SNPs} -b ${VCF_INVARIANT} > ${OUT_DIR}/temp_${OUT_PREFIX}.10x_snp_x_invariant.txt
bedtools intersect -sorted -wa -a ${VCF_10X_INDELs} -b ${VCF_INVARIANT} > ${OUT_DIR}/temp_${OUT_PREFIX}.10x_indel_x_invariant.txt
