#!/bin/bash

set -e
set -o pipefail

# Filter VCFs for morex-sample2 Nanopore data for the purposes of getting a set of diffs from reference

# Dependencies
module load bcftools/1.10.2
module load htslib/1.9
module load bedtools/2.29.2

# User provided input arguments
# VCF files
VCF="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3/Morex_ont_partsRefv3_90_renamed.vcf"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3/filtered"
# Output file prefix
OUT_PREFIX="morex_ont"
# Minimum number of reads supporting the variant call
MIN_SUPPORT="5"

# Uncallable regions
# This includes:
#   Regions where REF has stretches of N's
#   Repeat annotations
#   High copy regions (e.g., chloroplasts, mitochondria, rDNA repeats, centromere repeats, etc.)
#   Low complexity regions generated from JGI's BBMask, but includes 10% around gene models
UNCALLABLE="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/uncallable_regions/morex_v3_combined_uncallable.low_complexity.nochrUn.bed"

# Sniffles VCF to BED script
SNIFFLES_to_BED="/panfs/jay/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/01_snp_filtering/vcf_sniffles_long_read_to_bed.py"

#----------------
# Make output directories
mkdir -p ${OUT_DIR} ${OUT_DIR}/intermediates

# Remove sites that are homozygous reference, we want to identify sites that are differences from reference
# Exclude chrUn
# Remove sites below the min support threshold
bcftools filter --targets "^chrUn" -e "GT[*]='RR' | INFO/SUPPORT<${MIN_SUPPORT}" ${VCF} -O z -o ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.vcf.gz
# Index vcf
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.vcf.gz

# Separate indels from INV, DUP, and BND variants
# BND only
bcftools view -i 'INFO/SVTYPE="BND"' ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.vcf.gz -O v -o ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.BND.vcf
# DUP only
bcftools view -i 'INFO/SVTYPE="DUP"' ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.vcf.gz -O v -o ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.DUP.vcf
# INV only
bcftools view -i 'INFO/SVTYPE="INV"' ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.vcf.gz -O v -o ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.INV.vcf
# INS and DEL
bcftools view -i 'INFO/SVTYPE="INS" | INFO/SVTYPE="DEL"' ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.vcf.gz -O v -o ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.INDELs.vcf

# Remove SVs that overlap with uncallable regions
bedtools intersect -wa -v -header -a ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.INDELs.vcf -b ${UNCALLABLE} > ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.callable.INDELs.vcf
# Separate INS and DEL for exploration
bcftools view -i 'INFO/SVTYPE="INS"' ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.callable.INDELs.vcf -O v -o ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.callable.INS.vcf
bcftools view -i 'INFO/SVTYPE="DEL"' ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.callable.INDELs.vcf -O v -o ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.callable.DEL.vcf

# Convert Sniffles VCF to BED
python ${SNIFFLES_to_BED} ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.callable.INS.vcf > ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.callable.INS.bed
python ${SNIFFLES_to_BED} ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.callable.DEL.vcf > ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.callable.DEL.bed

# Bgzip and index (required for running samplot)
bgzip -c ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.callable.INS.bed > ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.callable.INS.bed.gz
tabix -p bed ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.callable.INS.bed.gz

bgzip -c ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.callable.DEL.bed > ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.callable.DEL.bed.gz
tabix -p bed ${OUT_DIR}/${OUT_PREFIX}.noHomRef.geSup${MIN_SUPPORT}.callable.DEL.bed.gz
