#!/bin/bash

set -e
set -o pipefail

# Dependencies
module load bcftools/1.10.2
module load htslib/1.9
module load bedtools/2.29.2

# User provided input arguments
VCF_M01="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/cutesv_calls/M01_ont.vcf"
VCF_M20="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/cutesv_calls/M20_ont.vcf"
VCF_M29="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/cutesv_calls/M29_ont.vcf"

OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/cutesv_calls/filtered"
OUT_PREFIX="mut_ont_cutesv"

SAMP1="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/cutesv_calls/filtered/sample_name-M01.txt"
SAMP2="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/cutesv_calls/filtered/sample_name-M20.txt"
SAMP3="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/cutesv_calls/filtered/sample_name-M29.txt"

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

#-----------------
mkdir -p ${OUT_DIR} ${OUT_DIR}/Intermediates

# Prepare file basenames
m01_bn=$(basename ${VCF_M01} .vcf)
m20_bn=$(basename ${VCF_M20} .vcf)
m29_bn=$(basename ${VCF_M29} .vcf)

# Remove chrUn and BND/DUP variants
# Remove homozygous reference sites
# Separate other variant types into separate files
# Include only sites that have at least one alt GT
# M01
# INS and DEL
bcftools view --targets "^chrUn" -e 'SVTYPE="BND" | SVTYPE="BND"' ${VCF_M01} | bcftools reheader --samples ${SAMP1} | bcftools view -e 'GT[*]="ref"' | bcftools view -i 'SVTYPE="INS" | SVTYPE="DEL"' -O z -o ${OUT_DIR}/Intermediates/${m01_bn}.INDELs.vcf.gz
# INV
bcftools view --targets "^chrUn" -e 'SVTYPE="BND" | SVTYPE="BND"' ${VCF_M01} | bcftools reheader --samples ${SAMP1} | bcftools view -e 'GT[*]="ref"' | bcftools view -i 'SVTYPE="INV"' -O z -o ${OUT_DIR}/Intermediates/${m01_bn}.INV.vcf.gz

# M20
# INS and DEL
bcftools view --targets "^chrUn" -e 'SVTYPE="BND" | SVTYPE="BND"' ${VCF_M20} | bcftools reheader --samples ${SAMP2} | bcftools view -e 'GT[*]="ref"' | bcftools view -i 'SVTYPE="INS" | SVTYPE="DEL"' -O z -o ${OUT_DIR}/Intermediates/${m20_bn}.INDELs.vcf.gz
# INV
bcftools view --targets "^chrUn" -e 'SVTYPE="BND" | SVTYPE="BND"' ${VCF_M20} | bcftools reheader --samples ${SAMP2} | bcftools view -e 'GT[*]="ref"' | bcftools view -i 'SVTYPE="INV"' -O z -o ${OUT_DIR}/Intermediates/${m20_bn}.INV.vcf.gz

# M29
# INS and DEL
bcftools view --targets "^chrUn" -e 'SVTYPE="BND" | SVTYPE="BND"' ${VCF_M29} | bcftools reheader --samples ${SAMP3} | bcftools view -e 'GT[*]="ref"' | bcftools view -i 'SVTYPE="INS" | SVTYPE="DEL"' -O z -o ${OUT_DIR}/Intermediates/${m29_bn}.INDELs.vcf.gz
# INV
bcftools view --targets "^chrUn" -e 'SVTYPE="BND" | SVTYPE="BND"' ${VCF_M29} | bcftools reheader --samples ${SAMP3} | bcftools view -e 'GT[*]="ref"' | bcftools view -i 'SVTYPE="INV"' -O z -o ${OUT_DIR}/Intermediates/${m29_bn}.INV.vcf.gz

# Index vcf
tabix -p vcf ${OUT_DIR}/Intermediates/${m01_bn}.INDELs.vcf.gz
tabix -p vcf ${OUT_DIR}/Intermediates/${m01_bn}.INV.vcf.gz
tabix -p vcf ${OUT_DIR}/Intermediates/${m20_bn}.INDELs.vcf.gz
tabix -p vcf ${OUT_DIR}/Intermediates/${m20_bn}.INV.vcf.gz
tabix -p vcf ${OUT_DIR}/Intermediates/${m29_bn}.INDELs.vcf.gz
tabix -p vcf ${OUT_DIR}/Intermediates/${m29_bn}.INV.vcf.gz

# Combine into single file
# INDEL
bcftools merge ${OUT_DIR}/Intermediates/${m01_bn}.INDELs.vcf.gz ${OUT_DIR}/Intermediates/${m20_bn}.INDELs.vcf.gz ${OUT_DIR}/Intermediates/${m29_bn}.INDELs.vcf.gz -O z -o ${OUT_DIR}/Intermediates/${OUT_PREFIX}.INDELs.vcf.gz
# INV
bcftools merge ${OUT_DIR}/Intermediates/${m01_bn}.INV.vcf.gz ${OUT_DIR}/Intermediates/${m20_bn}.INV.vcf.gz ${OUT_DIR}/Intermediates/${m29_bn}.INV.vcf.gz -O z -o ${OUT_DIR}/Intermediates/${OUT_PREFIX}.INV.vcf.gz

# Only include INDELs private to each line
# INDEL
bcftools view -i "COUNT(GT='alt')=1" ${OUT_DIR}/Intermediates/${OUT_PREFIX}.INDELs.vcf.gz -O z -o ${OUT_DIR}/${OUT_PREFIX}.private.INDELs.vcf.gz
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}.private.INDELs.vcf.gz
# INV
bcftools view -i "COUNT(GT='alt')=1" ${OUT_DIR}/Intermediates/${OUT_PREFIX}.INV.vcf.gz -O z -o ${OUT_DIR}/${OUT_PREFIX}.private.INV.vcf.gz

# Remove SVs that overlap with uncallable regions
# INDEL
bedtools intersect -wa -v -header -a ${OUT_DIR}/${OUT_PREFIX}.private.INDELs.vcf.gz -b ${UNCALLABLE} > ${OUT_DIR}/${OUT_PREFIX}.private.callable.INDELs.vcf
# INV
# Note: no output file was written since 0 INV remain after the uncallable regions filter
bedtools intersect -wa -v -a ${OUT_DIR}/${OUT_PREFIX}.private.INV.vcf.gz -b ${UNCALLABLE} | wc -l

# Remove SVs that overlap with diffs in ref of our Morex sample and Mascher et al. Morex
bedtools intersect -wa -v -header -a ${OUT_DIR}/${OUT_PREFIX}.private.callable.INDELs.vcf -b ${REF_DIFFS_10x_del} ${REF_DIFFS_ONT_DEL} ${REF_DIFFS_ONT_INS} ${REF_DIFFS_85xONT_DEL} ${REF_DIFFS_85xONT_INS} ${REF_DIFFS_PacBio_DEL} ${REF_DIFFS_PacBio_INS} | uniq > ${OUT_DIR}/${OUT_PREFIX}.private.callable.noRefDiffs.INDELs.vcf

# Separate INS and DEL
# Add tag "BasesToClosestVariant" to remove (mostly) consecutive variants since they are likely not de novo
bcftools view -i 'SVTYPE="INS"' ${OUT_DIR}/${OUT_PREFIX}.private.callable.noRefDiffs.INDELs.vcf | vcfdistance | bcftools view -e 'BasesToClosestVariant <= 10 & SVLEN > 10' | bgzip > ${OUT_DIR}/${OUT_PREFIX}.private.callable.noRefDiffs.final.INS.vcf.gz
bcftools view -i 'SVTYPE="DEL"' ${OUT_DIR}/${OUT_PREFIX}.private.callable.noRefDiffs.INDELs.vcf | vcfdistance | bcftools view -e 'BasesToClosestVariant <= 10 & SVLEN < -10' | bgzip > ${OUT_DIR}/${OUT_PREFIX}.private.callable.noRefDiffs.final.DEL.vcf.gz
# Index vcf files
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}.private.callable.noRefDiffs.final.INS.vcf.gz
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}.private.callable.noRefDiffs.final.DEL.vcf.gz

# Concatenate INS and DEL to final INDELs file
bcftools concat -a ${OUT_DIR}/${OUT_PREFIX}.private.callable.noRefDiffs.final.INS.vcf.gz ${OUT_DIR}/${OUT_PREFIX}.private.callable.noRefDiffs.final.DEL.vcf.gz -O v -o ${OUT_DIR}/${OUT_PREFIX}.private.callable.noRefDiffs.final.INDELs.vcf
