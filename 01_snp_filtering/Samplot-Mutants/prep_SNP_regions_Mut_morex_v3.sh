#!/bin/bash

set -e
set -o pipefail

# Prepare regions files for each variant type prior to running samplot.
# This script serves as a log of commands run to prepare files.
# Variant types we'll prepare are:
#   DEL, INS, INV, BND, DUP

# Dependencies
module load bedops_ML/2.4.38
module load bedtools/2.29.2

# User provided input arguments
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration_morex_v3"
# VCF files, one per variable
VCF_10x_M20="/scratch.global/liux1299/M20-2-2_phased_variants_snps_filtered_morex_v3.vcf"
VCF_10x_M29="/scratch.global/liux1299/M29-2-2_phased_variants_snps_filtered_morex_v3.vcf"
# Output file names
SAMPLE_NAME_M20="M20_10xGenomics_snps_morex_v3"
SAMPLE_NAME_M29="M29_10xGenomics_snps_morex_v3"
SNPS_OUTNAME="SNPs_M20_M29_10xGenomics_morex_v3"
# This bed file lists all intervals with 'Ns' in the reference genome
REF_Ns_BED="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/stretches_of_Ns/Barley_MorexV3_pseudomolecules_parts_missing.bed.gz"

#----------------
# Make output directory
mkdir -p ${OUT_DIR}

# Convert to BED format
# M20
vcf2bed < ${VCF_10x_M20} | cut -f 1,2,3 > ${OUT_DIR}/${SAMPLE_NAME_M20}.bed
# M29
vcf2bed < ${VCF_10x_M29} | cut -f 1,2,3 > ${OUT_DIR}/${SAMPLE_NAME_M29}.bed

# Remove stretches of Ns in REF
# also exclude chrUn
# M20
bedtools intersect -a ${OUT_DIR}/${SAMPLE_NAME_M20}.bed -b ${REF_Ns_BED} -v | grep -v "chrUn" > ${OUT_DIR}/SNPs_${SAMPLE_NAME_M20}_noMissing.bed
# M29
bedtools intersect -a ${OUT_DIR}/${SAMPLE_NAME_M29}.bed -b ${REF_Ns_BED} -v | grep -v "chrUn" > ${OUT_DIR}/SNPs_${SAMPLE_NAME_M29}_noMissing.bed

# Combine two BED files into a single file
cat ${OUT_DIR}/${SAMPLE_NAME_M20}.bed ${OUT_DIR}/${SAMPLE_NAME_M29}.bed | sort -k1,1 -k2,2n > ${OUT_DIR}/${SNPS_OUTNAME}.bed
bedtools merge -i ${OUT_DIR}/${SNPS_OUTNAME}.bed > ${OUT_DIR}/${SNPS_OUTNAME}_merged.bed
# Cleanup intermediate files
rm ${OUT_DIR}/${SNPS_OUTNAME}.bed
