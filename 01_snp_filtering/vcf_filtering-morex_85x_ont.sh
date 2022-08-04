#!/bin/bash

set -e
set -o pipefail

# Filter VCFs for Morex Nanopore data for the purposes of getting a set of diffs from reference

# Dependencies
module load bcftools/1.10.2
module load htslib/1.9
module load bedtools/2.29.2

# User provided input arguments
# VCF files
VCF="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_85x_ont/Morex_85X_sorted_renamed_parts.vcf"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_85x_ont/filtered"
# Output file prefix
OUT_PREFIX="morex_85x_ont"
# List of regions where REF has stretches of N's
REF_Ns_BED="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/stretches_of_Ns/Barley_MorexV3_pseudomolecules_parts_missing.bed"
# Minimum number of reads supporting the variant call
MIN_SUPPORT="3"
REPEAT_ANN="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/annotation/HvulgareMorex_702_V3.repeatmasked_assembly_V3.parts.gff3"

# High confidence genes
#HIGH_CONF_GENES="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/gene_annotation/Barley_Morex_V2_gene_annotation_PGSB.HC.parts.sorted.uniq.bed"

#----------------
# Make output directories
mkdir -p ${OUT_DIR} ${OUT_DIR}/intermediates

# Remove sites that are homozygous reference, we want to identify sites that are differences from reference
# We'll also exclude chrUn here too and AP017301.1 (since the sequences we are comparing to don't have AP017301.1)
bcftools view -e 'GT[*]="RR"' ${VCF} | grep -v "chrUn" | grep -v "AP017301.1" > ${OUT_DIR}/${OUT_PREFIX}_noHomRef.vcf

# Remove sites below the min support threshold
bcftools view -i "INFO/SUPPORT >= ${MIN_SUPPORT}" ${OUT_DIR}/${OUT_PREFIX}_noHomRef.vcf -O z -o ${OUT_DIR}/${OUT_PREFIX}_noHomRef_geSup${MIN_SUPPORT}.vcf.gz
# Index vcf
bcftools index ${OUT_DIR}/${OUT_PREFIX}_noHomRef_geSup${MIN_SUPPORT}.vcf.gz

# Remove SVs that overlap with repeat annotated regions and that overlap with stretches of Ns
bedtools intersect -wa -v -header -a ${OUT_DIR}/${OUT_PREFIX}_noHomRef_geSup${MIN_SUPPORT}.vcf.gz -b ${REPEAT_ANN} ${REF_Ns_BED} > ${OUT_DIR}/${OUT_PREFIX}_noHomRef_geSup${MIN_SUPPORT}.noRepeatOverlap.noRefNs.vcf
