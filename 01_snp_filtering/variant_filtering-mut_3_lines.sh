#!/bin/bash

set -e
set -o pipefail

# Filter 10x Genomics VCF file containing 3 mutated lines (includes SNPs, dels, and SVs)
# Output VCF file containing singletons

# Dependencies
module load bcftools/1.9
module load bedtools/2.29.2
module load vcftools_ML/0.1.16
# Required for vcf-annotate that is part of vcftools
export PERL5LIB=$PERL5LIB:/panfs/roc/groups/9/morrellp/public/Software/vcftools_ML-0.1.16/share/perl5

# User provided input arguments
RAW_VCF=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated/mut_3_lines_sorted.vcf
# Full filepath to output directory
OUT_DIR=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated/Filtered
# Output file prefix
PREFIX=mut_3_lines
# BED file containing sites that differ between 10x Morex and Morex reference
BED_EXCLUSION_LIST=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/morex-sample2/Filtered/morex-sample2_diffs_from_ref.bed
VCFTOOLS_CUSTOM_FILTER=~/GitHub/Barley_Mutated/01_snp_filtering/filters.txt

#----------------------------------------
# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR}

# Pass 1 filtering
# Filter out sites using 10x Genomics custom filters
bcftools filter -e 'FILTER=="10X_QUAL_FILTER" || FILTER=="10X_ALLELE_FRACTION_FILTER" || FILTER=="10X_PHASING_INCONSISTENT" || FILTER=="10X_HOMOPOLYMER_UNPHASED_INSERTION" || FILTER=="10X_RESCUED_MOLECULE_HIGH_DIVERSITY" || FILTER=="LOWQ"' ${RAW_VCF} > ${OUT_DIR}/${PREFIX}_filtered_pass1.vcf

# Pass 2
# Filter out sites with DP < 5
bcftools filter -e 'INFO/DP < 5' ${OUT_DIR}/${PREFIX}_filtered_pass1.vcf > ${OUT_DIR}/${PREFIX}_filtered_pass2.vcf

# Exclude sites that differ between 10x Morex and Morex reference
# Prepare header line first since bedtools doesn't preserve the header line
grep "#" ${OUT_DIR}/${PREFIX}_filtered_pass2.vcf > ${OUT_DIR}/${PREFIX}_filtered_no_morex_diffs.vcf
bedtools intersect -v -a ${OUT_DIR}/${PREFIX}_filtered_pass2.vcf -b ${BED_EXCLUSION_LIST} >> ${OUT_DIR}/${PREFIX}_filtered_no_morex_diffs.vcf

# Exclude non-unique variants (i.e., variants that are present in more than 1 of the mutated lines)
cat ${OUT_DIR}/${PREFIX}_filtered_no_morex_diffs.vcf | vcf-annotate -f ${VCFTOOLS_CUSTOM_FILTER} > ${OUT_DIR}/${PREFIX}_filtered_singleton_ann.vcf
# Create a file containing only singletons
grep "#" ${OUT_DIR}/${PREFIX}_filtered_singleton_ann.vcf > ${OUT_DIR}/${PREFIX}_filtered_singletons_only.vcf
grep -v "#" ${OUT_DIR}/${PREFIX}_filtered_singleton_ann.vcf | grep "SINGLETON" >> ${OUT_DIR}/${PREFIX}_filtered_singletons_only.vcf

# Pull out homozygous sites only
bcftools view -i 'GT[*]="hom"' ${OUT_DIR}/${PREFIX}_filtered_singletons_only.vcf > ${OUT_DIR}/${PREFIX}_filtered_hom_singletons_only.vcf
