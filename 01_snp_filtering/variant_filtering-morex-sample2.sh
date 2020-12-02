#!/bin/bash

set -e
set -o pipefail

# Filter 10x Genomics VCF file (includes SNPs, dels, and SVs) and identify differences between 10x Morex and Morex reference
#   This script outputs a VCF with differences from reference

# Dependencies
module load bcftools/1.9
module load bedops_ML/2.4.38

# User provided input arguments
RAW_VCF=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/morex-sample2/morex-sample2_concat_sorted.vcf
OUT_DIR=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/morex-sample2/Filtered
PREFIX=morex-sample2

# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR}

# First pass filtering
# Filter out sites using 10x Genomics custom filters
bcftools filter -e 'FILTER=="10X_QUAL_FILTER" || FILTER=="10X_ALLELE_FRACTION_FILTER" || FILTER=="10X_PHASING_INCONSISTENT" || FILTER=="10X_HOMOPOLYMER_UNPHASED_INSERTION" || FILTER=="10X_RESCUED_MOLECULE_HIGH_DIVERSITY" || FILTER=="LOWQ"' ${RAW_VCF} > ${OUT_DIR}/${PREFIX}_filtered_pass1.vcf

# Second pass filtering
# Filter out sites with DP < 5
bcftools filter -e 'INFO/DP < 5' ${OUT_DIR}/${PREFIX}_filtered_pass1.vcf > ${OUT_DIR}/${PREFIX}_filtered_pass2.vcf

# Third pass filtering
# Remove sites that are homozygous reference, we want to identify sites that are differences from reference
bcftools view -e 'GT[*]="RR"' ${OUT_DIR}/${PREFIX}_filtered_pass2.vcf > ${OUT_DIR}/${PREFIX}_diffs_from_ref.vcf

# Convert diffs from ref VCF to BED using bedops tool
# This file will be used as an exclusion list
vcf2bed < ${OUT_DIR}/${PREFIX}_diffs_from_ref.vcf | cut -f 1,2,3 > ${OUT_DIR}/${PREFIX}_diffs_from_ref.bed
