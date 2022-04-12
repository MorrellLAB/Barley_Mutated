#!/bin/bash

set -e
set -o pipefail

# Prepare regions files for each variant type prior to running samplot.
# This script serves as a log of commands run to prepare files.
# Variant types we'll prepare are:
#   DEL, INS, INV, BND, DUP

# Dependencies
module load bcftools/1.9
module load bedtools/2.29.2
module load bedops_ML/2.4.38
module load python3/3.8.3_anaconda2020.07_mamba

# User provided input arguments
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration_morex_v3"
# VCF files, one per variable
VCF_ONT="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3/filtered/morex_ont_noHomRef.vcf"
VCF_PacBio="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/pacbio_morex/pacbio_morex_v3/filtered/morex_pacbio_noHomRef.vcf"
SAMPLE_NAME="morex"
# Sample IDs, one per variable
# Must be in the same order as the VCF files but can be named anything and don't need
#   to match actual sample names in VCF files.
SAMPLE_ID_ONT="morex_ont"
SAMPLE_ID_PacBio="morex_pacbio"
# This bed file lists all intervals with 'Ns' in the reference genome
REF_Ns_BED="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/stretches_of_Ns/Barley_MorexV3_pseudomolecules_parts_missing.bed"

# Path to VCF to BED scripts
VCF_10x_TO_BED_SCRIPT="/panfs/roc/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/01_snp_filtering/vcf_10x_genomics_to_bed.py"
VCF_SNIFFLES_TO_BED_SCRIPT="/panfs/roc/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/01_snp_filtering/vcf_sniffles_long_read_to_bed.py"

#----------------
# Make output directory
mkdir -p ${OUT_DIR}

# DEL: pull variant type from VCFs
# Nanopore VCF, called by Sniffles2
bcftools view -i 'INFO/SVTYPE=="DEL"' ${VCF_ONT} > ${OUT_DIR}/${SAMPLE_ID_ONT}_DEL.vcf
# PacBio VCF, called by Sniffles2
bcftools view -i 'INFO/SVTYPE=="DEL"' ${VCF_PacBio} > ${OUT_DIR}/${SAMPLE_ID_PacBio}_DEL.vcf
# Convert to BED format
python3 ${VCF_SNIFFLES_TO_BED_SCRIPT} ${OUT_DIR}/${SAMPLE_ID_ONT}_DEL.vcf > ${OUT_DIR}/${SAMPLE_ID_ONT}_DEL.bed
python3 ${VCF_SNIFFLES_TO_BED_SCRIPT} ${OUT_DIR}/${SAMPLE_ID_PacBio}_DEL.vcf > ${OUT_DIR}/${SAMPLE_ID_PacBio}_DEL.bed
# Concat BED files into one
cat ${OUT_DIR}/${SAMPLE_ID_ONT}_DEL.bed ${OUT_DIR}/${SAMPLE_ID_PacBio}_DEL.bed | sort -k1,1 -k2,2n > ${OUT_DIR}/DEL_${SAMPLE_NAME}_ONT_PacBio_sorted.bed
# Remove stretches of Ns in REF
bedtools intersect -a ${OUT_DIR}/DEL_${SAMPLE_NAME}_ONT_PacBio_sorted.bed -b ${REF_Ns_BED} -v > ${OUT_DIR}/DEL_${SAMPLE_NAME}_ONT_PacBio_noMissing.bed
# Merge BED files, exclude chrUn
#bedtools merge -i ${OUT_DIR}/DEL_${SAMPLE_NAME}_10x_ONT_PacBio_noMissing.bed | grep -v "chrUn" > ${OUT_DIR}/DEL_${SAMPLE_NAME}_10x_ONT_PacBio_merged.bed

# INS: pull variant type from VCFs
# Nanopore VCF, called by Sniffles2
bcftools view -i 'INFO/SVTYPE=="INS"' ${VCF_ONT} > ${OUT_DIR}/${SAMPLE_ID_ONT}_INS.vcf
# PacBio VCF, called by Sniffles2
bcftools view -i 'INFO/SVTYPE=="INS"' ${VCF_PacBio} > ${OUT_DIR}/${SAMPLE_ID_PacBio}_INS.vcf
# Convert to BED format
python3 ${VCF_SNIFFLES_TO_BED_SCRIPT} ${OUT_DIR}/${SAMPLE_ID_ONT}_INS.vcf > ${OUT_DIR}/${SAMPLE_ID_ONT}_INS.bed
python3 ${VCF_SNIFFLES_TO_BED_SCRIPT} ${OUT_DIR}/${SAMPLE_ID_PacBio}_INS.vcf > ${OUT_DIR}/${SAMPLE_ID_PacBio}_INS.bed
# Concat BED files into one
cat ${OUT_DIR}/${SAMPLE_ID_ONT}_INS.bed ${OUT_DIR}/${SAMPLE_ID_PacBio}_INS.bed | sort -k1,1 -k2,2n > ${OUT_DIR}/INS_${SAMPLE_NAME}_ONT_PacBio_sorted.bed
# Remove stretches of Ns in REF
bedtools intersect -a ${OUT_DIR}/INS_${SAMPLE_NAME}_ONT_PacBio_sorted.bed -b ${REF_Ns_BED} -v > ${OUT_DIR}/INS_${SAMPLE_NAME}_ONT_PacBio_noMissing.bed
# Merge BED files, exclude chrUn
#bedtools merge -i ${OUT_DIR}/INS_${SAMPLE_NAME}_10x_ONT_PacBio_noMissing.bed | grep -v "chrUn" > ${OUT_DIR}/INS_${SAMPLE_NAME}_10x_ONT_PacBio_merged.bed

# INV: pull variant type from VCFs
# Nanopore VCF, called by Sniffles2
bcftools view -i 'INFO/SVTYPE=="INV"' ${VCF_ONT} > ${OUT_DIR}/${SAMPLE_ID_ONT}_INV.vcf
# PacBio VCF, called by Sniffles2
bcftools view -i 'INFO/SVTYPE=="INV"' ${VCF_PacBio} > ${OUT_DIR}/${SAMPLE_ID_PacBio}_INV.vcf
# Convert to BED format
python3 ${VCF_SNIFFLES_TO_BED_SCRIPT} ${OUT_DIR}/${SAMPLE_ID_ONT}_INV.vcf > ${OUT_DIR}/${SAMPLE_ID_ONT}_INV.bed
python3 ${VCF_SNIFFLES_TO_BED_SCRIPT} ${OUT_DIR}/${SAMPLE_ID_PacBio}_INV.vcf > ${OUT_DIR}/${SAMPLE_ID_PacBio}_INV.bed
# Concat BED files into one
cat ${OUT_DIR}/${SAMPLE_ID_ONT}_INV.bed ${OUT_DIR}/${SAMPLE_ID_PacBio}_INV.bed | sort -k1,1 -k2,2n > ${OUT_DIR}/INV_${SAMPLE_NAME}_ONT_PacBio_sorted.bed
