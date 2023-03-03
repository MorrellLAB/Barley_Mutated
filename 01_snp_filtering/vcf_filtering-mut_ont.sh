#!/bin/bash

set -e
set -o pipefail

# Dependencies
module load bcftools/1.10.2
module load htslib/1.9
module load bedtools/2.29.2
module load python3/3.8.3_anaconda2020.07_mamba

# User provided input arguments
VCF="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/mut_ont.vcf.gz"
# Output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/filtered"
# Output file prefix
OUT_PREFIX="mut_ont"
# Minimum number of reads supporting the variant call
MIN_SUPPORT="3"

# Uncallable regions
# This includes:
#   Regions where REF has stretches of N's
#   Repeat annotations
#   High copy regions (e.g., chloroplasts, mitochondria, rDNA repeats, centromere repeats, etc.)
#   Low complexity regions generated from JGI's BBMask, but includes 10% around gene models
UNCALLABLE="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/uncallable_regions/morex_v3_combined_uncallable.low_complexity.nochrUn.bed"

# Differences from Morex reference
REF_DIFFS_10x_del="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/diffs_from_ref/morex-sample2_dels_diffs_from_ref.noBND.bed"
REF_DIFFS_ONT_DEL="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3/filtered/morex_ont.noHomRef.geSup5.callable.DEL.bed"
REF_DIFFS_ONT_INS="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3/filtered/morex_ont.noHomRef.geSup5.callable.INS.bed"
REF_DIFFS_85xONT_DEL="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_85x_ont/filtered/parts_pos/morex_85x_ont.noHomRef.geSup5.callable.parts.DEL.bed"
REF_DIFFS_85xONT_INS="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_85x_ont/filtered/parts_pos/morex_85x_ont.noHomRef.geSup5.callable.parts.INS.bed"
REF_DIFFS_PacBio_DEL="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/pacbio_morex/pacbio_morex_v3/filtered/morex_pacbio.noHomRef.geSup5.callable.DEL.bed"
REF_DIFFS_PacBio_INS="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/pacbio_morex/pacbio_morex_v3/filtered/morex_pacbio.noHomRef.geSup5.callable.INS.bed"

# Sniffles VCF to BED script
SNIFFLES_to_BED="/panfs/jay/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/01_snp_filtering/vcf_sniffles_long_read_to_bed.py"

#----------------
# Make output directories
mkdir -p ${OUT_DIR} ${OUT_DIR}/intermediates

# Include only sites that have at least one alt GT
# and SVs private to each line
# Filter on minimum support
bcftools view -i "GT[*]='alt' & COUNT(GT='alt')=1" ${VCF} | bcftools view -i "INFO/SUPPORT>=${MIN_SUPPORT}" -O z -o ${OUT_DIR}/${OUT_PREFIX}.private.geSup${MIN_SUPPORT}.vcf.gz
# Index vcf
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}.private.geSup${MIN_SUPPORT}.vcf.gz

# Separate indels from INV, DUP, and BND variants
# BND only
bcftools view -i 'INFO/SVTYPE="BND"' ${OUT_DIR}/${OUT_PREFIX}.private.geSup${MIN_SUPPORT}.vcf.gz -O v -o ${OUT_DIR}/${OUT_PREFIX}.private.geSup${MIN_SUPPORT}.BND.vcf
# DUP only
bcftools view -i 'INFO/SVTYPE="DUP"' ${OUT_DIR}/${OUT_PREFIX}.private.geSup${MIN_SUPPORT}.vcf.gz -O v -o ${OUT_DIR}/${OUT_PREFIX}.private.geSup${MIN_SUPPORT}.DUP.vcf
# INV only
bcftools view -i 'INFO/SVTYPE="INV"' ${OUT_DIR}/${OUT_PREFIX}.private.geSup${MIN_SUPPORT}.vcf.gz -O v -o ${OUT_DIR}/${OUT_PREFIX}.private.geSup${MIN_SUPPORT}.INV.vcf
# INS and DEL
bcftools view -i 'INFO/SVTYPE="INS" | INFO/SVTYPE="DEL"' ${OUT_DIR}/${OUT_PREFIX}.private.geSup${MIN_SUPPORT}.vcf.gz -O v -o ${OUT_DIR}/${OUT_PREFIX}.private.geSup${MIN_SUPPORT}.INDELs.vcf

# Remove SVs that overlap with uncallable regions
# For DUP and INV, zero SVs remain after considering uncallable regions
# So, we'll just proceed with indels
bedtools intersect -wa -v -header -a ${OUT_DIR}/${OUT_PREFIX}.private.geSup${MIN_SUPPORT}.INDELs.vcf -b ${UNCALLABLE} > ${OUT_DIR}/${OUT_PREFIX}.private.geSup${MIN_SUPPORT}.callable.INDELs.vcf
# Separate INS and DEL for exploration
bcftools view -i 'INFO/SVTYPE="INS"' ${OUT_DIR}/${OUT_PREFIX}.private.geSup${MIN_SUPPORT}.callable.INDELs.vcf -O v -o ${OUT_DIR}/${OUT_PREFIX}.private.geSup${MIN_SUPPORT}.callable.INS.vcf
bcftools view -i 'INFO/SVTYPE="DEL"' ${OUT_DIR}/${OUT_PREFIX}.private.geSup${MIN_SUPPORT}.callable.INDELs.vcf -O v -o ${OUT_DIR}/${OUT_PREFIX}.private.geSup${MIN_SUPPORT}.callable.DEL.vcf

# Remove SVs that overlap with diffs in ref of our Morex sample and Mascher et al. Morex
bedtools intersect -wa -v -header -a ${OUT_DIR}/${OUT_PREFIX}.private.geSup${MIN_SUPPORT}.callable.INDELs.vcf -b ${REF_DIFFS_10x_del} ${REF_DIFFS_ONT_DEL} ${REF_DIFFS_ONT_INS} ${REF_DIFFS_85xONT_DEL} ${REF_DIFFS_85xONT_INS} ${REF_DIFFS_PacBio_DEL} ${REF_DIFFS_PacBio_INS} > ${OUT_DIR}/${OUT_PREFIX}.private.geSup${MIN_SUPPORT}.callable.noRefDiffs.INDELs.vcf
# Separate INS and DEL for exploration
bcftools view -i 'INFO/SVTYPE="INS"' ${OUT_DIR}/${OUT_PREFIX}.private.geSup${MIN_SUPPORT}.callable.noRefDiffs.INDELs.vcf -O v -o ${OUT_DIR}/${OUT_PREFIX}.private.geSup${MIN_SUPPORT}.callable.noRefDiffs.INS.vcf
bcftools view -i 'INFO/SVTYPE="DEL"' ${OUT_DIR}/${OUT_PREFIX}.private.geSup${MIN_SUPPORT}.callable.noRefDiffs.INDELs.vcf -O v -o ${OUT_DIR}/${OUT_PREFIX}.private.geSup${MIN_SUPPORT}.callable.noRefDiffs.DEL.vcf

# Prepare BED format
${SNIFFLES_to_BED} ${OUT_DIR}/${OUT_PREFIX}.private.geSup${MIN_SUPPORT}.callable.noRefDiffs.INDELs.vcf > ${OUT_DIR}/${OUT_PREFIX}.private.geSup${MIN_SUPPORT}.callable.noRefDiffs.INDELs.bed
