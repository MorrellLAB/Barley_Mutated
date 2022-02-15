#!/bin/bash

set -e
set -o pipefail

# Filter 10x Genomics VCF file (includes SNPs, dels, and SVs) and identify differences between 10x Morex and Morex reference
#   This script outputs a VCF with differences from reference

# Dependencies
module load bcftools/1.9
module load bedops_ML/2.4.38
module load bedtools/2.29.2
module load python3/3.8.3_anaconda2020.07_mamba
# Export path to directory containing custom script for converting 10x Genomics VCF
#   containing dels, dups, and SVs where the end position is in the INFO field
export PATH=${PATH}:/panfs/roc/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/01_snp_filtering

# User provided input arguments
RAW_VCF="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/morex-sample2/morex-sample2_concat_sorted.vcf"
BEDPE="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/morex-sample2/morex-sample2_large_sv_calls.bedpe"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/morex-sample2/Filtered"
PREFIX="morex-sample2"

# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR}

# First pass filtering
# Filter out sites using 10x Genomics custom filters
# See here for info on filters:
#   https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/vcf
# Also remove sites where there are unplaced scaffolds (chrUn)
bcftools filter -e 'FILTER=="10X_QUAL_FILTER" || FILTER=="10X_ALLELE_FRACTION_FILTER" || FILTER=="10X_PHASING_INCONSISTENT" || FILTER=="10X_HOMOPOLYMER_UNPHASED_INSERTION" || FILTER=="10X_RESCUED_MOLECULE_HIGH_DIVERSITY" || FILTER=="LOWQ"' ${RAW_VCF} | grep -v "chrUn" > ${OUT_DIR}/${PREFIX}_filtered_pass1.vcf

# Second pass filtering
# A) Separate variants with the DP annotation from ones that don't have the DP annotation
# This is unique to the 10x Genomics VCF format where phased variants (i.e., SNPs, some indels)
#   have a DP annotation, but the dels and SVs VCFs don't have the DP annotation
#   so we have to filter them separately.
# B) For variants with DP annotation, keep sites where DP > 5 and DP < 78
# In this case there is only one sample so the INFO/DP and FORMAT DP values when they both exist
#   are the same.
# C) Remove sites that are homozygous reference, we want to identify sites that are differences from reference
bcftools view -i 'INFO/DP!="."' ${OUT_DIR}/${PREFIX}_filtered_pass1.vcf | bcftools filter -i 'INFO/DP > 5 && INFO/DP < 78' | bcftools view -e 'GT[*]="RR"' > ${OUT_DIR}/${PREFIX}_filtered.dp_ann.vcf
# Pull out variants without DP annotation that are NOT BND (SV breakend sets)
bcftools view -i 'INFO/DP=="."' ${OUT_DIR}/${PREFIX}_filtered_pass1.vcf | bcftools view -e 'INFO/SVTYPE=="BND" || INFO/SVTYPE2=="BND"' | bcftools view -e 'GT[*]="RR"' > ${OUT_DIR}/${PREFIX}_filtered.SVs_not_BND.vcf
# Pull out BND (SV breakend sets)
#   These need to be dealt with differently when converting to BED format
bcftools view -i 'INFO/DP=="." && (INFO/SVTYPE=="BND" | INFO/SVTYPE2=="BND")' ${OUT_DIR}/${PREFIX}_filtered_pass1.vcf | bcftools view -e 'GT[*]="RR"' > ${OUT_DIR}/${PREFIX}_filtered.SVs_BND.vcf

# Convert diffs from ref VCF to BED using bedops tool for typical VCF format (phased variants)
# and use a custom script for 10x Genomics specific dels, dups, and SVs formatting
# These files will be used as exclusion lists
vcf2bed < ${OUT_DIR}/${PREFIX}_filtered.dp_ann.vcf | cut -f 1,2,3 > ${OUT_DIR}/${PREFIX}_diffs_from_ref.dp_ann.bed
# Custom script for 10x Genomics specific format that are NOT BND (breakend sets)
vcf_10x_genomics_to_bed.py ${OUT_DIR}/${PREFIX}_filtered.SVs_not_BND.vcf > ${OUT_DIR}/${PREFIX}_diffs_from_ref.SVs_not_BND.bed
# BND sites VCF to BEDPE
bedtools intersect -wa -a ${BEDPE} -b ${OUT_DIR}/${PREFIX}_filtered.SVs_BND.vcf | sort -uV -k1,3 > ${OUT_DIR}/${PREFIX}_filtered.SVs_BND.bedpe
