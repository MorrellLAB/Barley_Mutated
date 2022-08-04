#!/bin/bash

set -e
set -o pipefail

# Filter anno.vcf.gz output from smoove called for all WGS short read data

# Dependencies
module load bcftools_PML/1.15.1
module load htslib/1.9
module load bedtools/2.29.2

# User provided input arguments
# VCF files
VCF="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/smoove_processing/results_smoove/results_genotyped/mut_barley_cohort_sites.smoove.square.anno.vcf.gz"
# Where do we output our files?
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/smoove_processing/results_smoove/results_genotyped/filtered"
# Output file prefix
PREFIX="mut_barley_cohort"

# List of regions where REF has stretches of N's
REF_Ns_BED="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/stretches_of_Ns/Barley_MorexV3_pseudomolecules_parts_missing.bed"
# Repeat annotations
REPEAT_ANN="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/annotation/HvulgareMorex_702_V3.repeatmasked_assembly_V3.parts.gff3"
# High copy regions (e.g., chloroplasts, mitochondria, rDNA repeats, centromere repeats, etc.)
HIGH_COPY_BED="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/high_copy_regions/Morex_v3_high_copy_uniq.parts.bed"

#---------------------
# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR}

# Separate BND variants from other variant types
vcf_dir=$(dirname ${VCF})
vcf_prefix=$(basename ${VCF} .vcf.gz)
# File of only BND variants
bcftools view -i 'INFO/SVTYPE="BND"' ${VCF} -O v -o ${vcf_dir}/${vcf_prefix}.BND_only.vcf
# File of other SVs, excluded BNDs
bcftools view -e 'INFO/SVTYPE="BND"' ${VCF} -O v -o ${vcf_dir}/${vcf_prefix}.noBND.vcf

# Filter on allelic balance (AB) proportion, expectation is it is close to 0.5.
#   Sets heterozygotes to missing if the AB deviates too much from 0.5.
#   Use '&' to select for sites where conditions are satisfied within the same sample
# Set per sample DP to missing if below threshold and remove chrUn
bcftools +setGT ${vcf_dir}/${vcf_prefix}.noBND.vcf -- -t q -n "." -i 'GT="het" & FMT/AB<0.3 | FMT/AB>0.7' | bcftools +setGT -- -t q -n "." -i 'FMT/DP<5' | bcftools filter --targets "^chrUn" -O z -o ${OUT_DIR}/${PREFIX}.filtAB-DP.vcf.gz
# Index vcf
tabix -p vcf ${OUT_DIR}/${PREFIX}.filtAB-DP.vcf.gz

# Remove SVs that overlap with repeat annotated regions and that overlap with stretches of Ns
bedtools intersect -wa -v -header -a ${OUT_DIR}/${PREFIX}.filtAB-DP.vcf.gz -b ${REPEAT_ANN} ${HIGH_COPY_BED} ${REF_Ns_BED} | bgzip > ${OUT_DIR}/${PREFIX}.filtAB-DP.noRepeatOverlap.noRefNs.vcf.gz
# Index vcf
tabix -p vcf ${OUT_DIR}/${PREFIX}.filtAB-DP.noRepeatOverlap.noRefNs.vcf.gz
