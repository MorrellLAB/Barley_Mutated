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
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration/samplot_prep"
# VCF files, one per variable
VCF_10x="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/morex-sample2/Filtered/morex-sample2_filtered_pass1.vcf"
VCF_ONT="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_nanopore_V2_partsRef/Morex_nanopore_V2_partsRef_90.vcf"
VCF_PacBio="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/pacbio_morex/pacbio_morex_v2/Morex_pacbio_90.vcf"
SAMPLE_NAME="morex"
# Sample IDs, one per variable
# Must be in the same order as the VCF files but can be named anything and don't need
#   to match actual sample names in VCF files.
SAMPLE_ID_10x="Morex_10x"
SAMPLE_ID_ONT="Morex_nanopore"
SAMPLE_ID_PacBio="Morex_pacbio"
# This bed file lists all intervals with 'Ns' in the reference genome
REF_Ns_BED="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration/Barley_Morex_V2_pseudomolecules_parts_missing.bed"

# Path to VCF to BED scripts
VCF_10x_TO_BED_SCRIPT="/panfs/roc/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/00_sequence_processing/vcf_10x_genomics_to_bed.py"
VCF_SNIFFLES_TO_BED_SCRIPT="/panfs/roc/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/00_sequence_processing/vcf_sniffles_long_read_to_bed.py"

#----------------
# Make output directory
mkdir -p ${OUT_DIR}

# DEL: pull variant type from VCFs
# 10x Genomics VCF
bcftools view -i 'INFO/SVTYPE=="DEL"' ${VCF_10x} > ${OUT_DIR}/${SAMPLE_ID_10x}_DEL.vcf
# Nanopore VCF, called by Sniffles2
bcftools view -i 'INFO/SVTYPE=="DEL"' ${VCF_ONT} > ${OUT_DIR}/${SAMPLE_ID_ONT}_DEL.vcf
# PacBio VCF, called by Sniffles2
bcftools view -i 'INFO/SVTYPE=="DEL"' ${VCF_PacBio} > ${OUT_DIR}/${SAMPLE_ID_PacBio}_DEL.vcf
# Convert to BED format
python3 ${VCF_10x_TO_BED_SCRIPT} ${OUT_DIR}/${SAMPLE_ID_10x}_DEL.vcf > ${OUT_DIR}/${SAMPLE_ID_10x}_DEL.bed
python3 ${VCF_SNIFFLES_TO_BED_SCRIPT} ${OUT_DIR}/${SAMPLE_ID_ONT}_DEL.vcf > ${OUT_DIR}/${SAMPLE_ID_ONT}_DEL.bed
python3 ${VCF_SNIFFLES_TO_BED_SCRIPT} ${OUT_DIR}/${SAMPLE_ID_PacBio}_DEL.vcf > ${OUT_DIR}/${SAMPLE_ID_PacBio}_DEL.bed
# Concat BED files into one
cat ${OUT_DIR}/${SAMPLE_ID_10x}_DEL.bed ${OUT_DIR}/${SAMPLE_ID_ONT}_DEL.bed ${OUT_DIR}/${SAMPLE_ID_PacBio}_DEL.bed | sort -k1,1 -k2,2n > ${OUT_DIR}/DEL_${SAMPLE_NAME}_10x_ONT_PacBio_sorted.bed
# Remove stretches of Ns in REF
bedtools intersect -a ${OUT_DIR}/DEL_${SAMPLE_NAME}_10x_ONT_PacBio_sorted.bed -b ${REF_Ns_BED} -v > ${OUT_DIR}/DEL_${SAMPLE_NAME}_10x_ONT_PacBio_noMissing.bed
# Merge BED files, exclude chrUn
bedtools merge -i ${OUT_DIR}/DEL_${SAMPLE_NAME}_10x_ONT_PacBio_noMissing.bed | grep -v "chrUn" > ${OUT_DIR}/DEL_${SAMPLE_NAME}_10x_ONT_PacBio_merged.bed

# INS: pull variant type from VCFs
# 10x Genomics VCF
# Note: bcftools can't work with bcftools view -i 'INFO/TYPE="ins"' ${VCF_10x}
#   because TYPE is a key word for bcftools. See this issue for details:
#   https://github.com/tgen/vcfMerger2/issues/18
# So, we'll prep this VCF a different way
grep "#" ${VCF_10x} > ${OUT_DIR}/${SAMPLE_ID_10x}_INS.vcf
grep -v "#" ${VCF_10x} | grep "TYPE=ins" >> ${OUT_DIR}/${SAMPLE_ID_10x}_INS.vcf
# Nanopore VCF, called by Sniffles2
bcftools view -i 'INFO/SVTYPE=="INS"' ${VCF_ONT} > ${OUT_DIR}/${SAMPLE_ID_ONT}_INS.vcf
# PacBio VCF, called by Sniffles2
bcftools view -i 'INFO/SVTYPE=="INS"' ${VCF_PacBio} > ${OUT_DIR}/${SAMPLE_ID_PacBio}_INS.vcf
# Convert to BED format
vcf2bed < ${OUT_DIR}/${SAMPLE_ID_10x}_INS.vcf | cut -f 1-3 > ${OUT_DIR}/${SAMPLE_ID_10x}_INS.bed
python3 ${VCF_SNIFFLES_TO_BED_SCRIPT} ${OUT_DIR}/${SAMPLE_ID_ONT}_INS.vcf > ${OUT_DIR}/${SAMPLE_ID_ONT}_INS.bed
python3 ${VCF_SNIFFLES_TO_BED_SCRIPT} ${OUT_DIR}/${SAMPLE_ID_PacBio}_INS.vcf > ${OUT_DIR}/${SAMPLE_ID_PacBio}_INS.bed
# Concat BED files into one
cat ${OUT_DIR}/${SAMPLE_ID_10x}_INS.bed ${OUT_DIR}/${SAMPLE_ID_ONT}_INS.bed ${OUT_DIR}/${SAMPLE_ID_PacBio}_INS.bed | sort -k1,1 -k2,2n > ${OUT_DIR}/INS_${SAMPLE_NAME}_10x_ONT_PacBio_sorted.bed
# Remove stretches of Ns in REF
bedtools intersect -a ${OUT_DIR}/INS_${SAMPLE_NAME}_10x_ONT_PacBio_sorted.bed -b ${REF_Ns_BED} -v > ${OUT_DIR}/INS_${SAMPLE_NAME}_10x_ONT_PacBio_noMissing.bed
# Merge BED files, exclude chrUn
bedtools merge -i ${OUT_DIR}/INS_${SAMPLE_NAME}_10x_ONT_PacBio_noMissing.bed | grep -v "chrUn" > ${OUT_DIR}/INS_${SAMPLE_NAME}_10x_ONT_PacBio_merged.bed

# INV: pull variant type from VCFs
# 10x Genomics VCF
bcftools view -i 'INFO/SVTYPE=="INV"' ${VCF_10x} > ${OUT_DIR}/${SAMPLE_ID_10x}_INV.vcf

