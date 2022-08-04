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
VCF="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3/Morex_ont_partsRefv3_90_renamed.vcf"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3/filtered"
# Output file prefix
OUT_PREFIX="morex_ont"
# List of regions where REF has stretches of N's
REF_Ns_BED="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/stretches_of_Ns/Barley_MorexV3_pseudomolecules_parts_missing.bed"
# Minimum number of reads supporting the variant call
MIN_SUPPORT="1"
REPEAT_ANN="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/annotation/HvulgareMorex_702_V3.repeatmasked_assembly_V3.parts.gff3"
# High copy regions (e.g., chloroplasts, mitochondria, rDNA repeats, centromere repeats, etc.)
HIGH_COPY_BED="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/high_copy_regions/Morex_v3_high_copy_uniq.parts.bed"
# High confidence genes
#HIGH_CONF_GENES="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/gene_annotation/Barley_Morex_V2_gene_annotation_PGSB.HC.parts.sorted.uniq.bed"

#----------------
# Make output directories
mkdir -p ${OUT_DIR} ${OUT_DIR}/intermediates

# Separate BND variants from other variant types
vcf_dir=$(dirname ${VCF})
vcf_prefix=$(basename ${VCF} .vcf)
# File of only BND variants
bcftools view -i 'INFO/SVTYPE="BND"' ${VCF} -O v -o ${vcf_dir}/${vcf_prefix}.BND_only.vcf
# File of other SVs, excluded BNDs
bcftools view -e 'INFO/SVTYPE="BND"' ${VCF} -O v -o ${vcf_dir}/${vcf_prefix}.noBND.vcf

# Remove sites that are homozygous reference, we want to identify sites that are differences from reference
# We'll also exclude chrUn here too
bcftools view -e 'GT[*]="RR"' ${vcf_dir}/${vcf_prefix}.noBND.vcf | grep -v "chrUn" > ${OUT_DIR}/${OUT_PREFIX}_noHomRef.vcf
# Compress
bgzip ${OUT_DIR}/${OUT_PREFIX}_noHomRef.vcf
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}_noHomRef.vcf.gz

# Remove sites below the min support threshold
bcftools view -i "INFO/SUPPORT >= ${MIN_SUPPORT}" ${OUT_DIR}/${OUT_PREFIX}_noHomRef.vcf.gz -O z -o ${OUT_DIR}/${OUT_PREFIX}_noHomRef_geSup${MIN_SUPPORT}.vcf.gz
# Index vcf
bcftools index ${OUT_DIR}/${OUT_PREFIX}_noHomRef_geSup${MIN_SUPPORT}.vcf.gz

# Remove SVs that overlap with repeat annotated regions and that overlap with stretches of Ns
bedtools intersect -wa -v -header -a ${OUT_DIR}/${OUT_PREFIX}_noHomRef_geSup${MIN_SUPPORT}.vcf.gz -b ${REPEAT_ANN} ${HIGH_COPY_BED} ${REF_Ns_BED} > ${OUT_DIR}/${OUT_PREFIX}_noHomRef_geSup${MIN_SUPPORT}.noRepeatOverlap.noRefNs.vcf
