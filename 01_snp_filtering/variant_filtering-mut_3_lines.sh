#!/bin/bash

set -e
set -o pipefail

# Filter 10x Genomics VCF file containing 3 mutated lines (includes SNPs, dels, and SVs)

# Dependencies
module load bcftools/1.9

# User provided input arguments
RAW_VCF=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated/mut_3_lines_sorted.vcf
OUT_DIR=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated/Filtered
PREFIX=mut_3_lines

# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR}

# First pass filtering
# Filter out sites using 10x Genomics custom filters
bcftools filter -e 'FILTER=="10X_QUAL_FILTER" || FILTER=="10X_ALLELE_FRACTION_FILTER" || FILTER=="10X_HOMOPOLYMER_UNPHASED_INSERTION" || FILTER=="10X_RESCUED_MOLECULE_HIGH_DIVERSITY" || FILTER=="LOWQ"' ${RAW_VCF} > ${OUT_DIR}/${PREFIX}_filtered_pass1.vcf

# Second pass
# Filter out sites with DP < 5
bcftools filter -e 'INFO/DP < 5' ${OUT_DIR}/${PREFIX}_filtered_pass1.vcf > ${OUT_DIR}/${PREFIX}_filtered_pass2.vcf

# Third pass
# Filter out sites with greater than 25% missing entries
bcftools filter -e 'F_PASS(GT="mis") > 0.25' ${OUT_DIR}/${PREFIX}_filtered_pass2.vcf > ${OUT_DIR}/${PREFIX}_filtered_pass3.vcf
