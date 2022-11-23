#!/bin/bash

# This script prepares uncallable regions for filtering.
# Some low complexity regions overlap with gene annotations. We want to keep
#   SVs called in gene models, so we'll remove the regions if there is overlap.

# Dependencies
module load bedtools/2.29.2

# User provided input arguments
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/uncallable_regions"
REF_FAI="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta.fai"
# BED file
GENE_ANN="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/gene_annotation/Hv_Morex.pgsb.Jul2020.sorted.parts.nochrUn.bed"
RESIZE_FRACTION="0.1"

# Combined BED files will have this output file prefix
OUT_PREFIX="morex_v3_combined_uncallable"

# Uncallable regions
# List of regions where REF has stretches of N's
REF_Ns_BED="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/stretches_of_Ns/Barley_MorexV3_pseudomolecules_parts_missing.bed"
# Repeat annotations
REPEAT_ANN="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/annotation/HvulgareMorex_702_V3.repeatmasked_assembly_V3.parts.bed"
# High copy regions (e.g., chloroplasts, mitochondria, rDNA repeats, centromere repeats, etc.)
HIGH_COPY_BED="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/high_copy_regions/Morex_v3_high_copy_uniq.parts.bed"
# Low complexity regions generated from JGI's BBMask
LOW_COMPLEXITY="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/entropy_masked/Barley_MorexV3_pseudomolecules_parts.entropy_0.7_masked.bed"

#-----------------
mkdir -p ${OUT_DIR}

# For exploration to see what intersects with gene annotations file
# Will help us decide which regions should be excluded
# bedtools intersect -wa -a ${REF_Ns_BED} -b ${GENE_ANN} | sort -uV -k1,1 -k2,2n
# bedtools intersect -wa -a ${REPEAT_ANN} -b ${GENE_ANN} | sort -uV -k1,1 -k2,2n
# bedtools intersect -wa -a ${HIGH_COPY_BED} -b ${GENE_ANN} | sort -uV -k1,1 -k2,2n
# bedtools intersect -wb -a ${HIGH_COPY_BED} -b ${GENE_ANN} | sort -uV -k1,1 -k2,2n
# bedtools intersect -wa -a ${LOW_COMPLEXITY} -b ${GENE_ANN} | sort -uV -k1,1 -k2,2n
# bedtools intersect -wb -a ${LOW_COMPLEXITY} -b ${GENE_ANN} | sort -uV -k1,1 -k2,2n

# Prepare output file prefix
ann_prefix=$(basename ${GENE_ANN} .bed)
# We'll preserve sequence 10% up or downstream of gene annotations
# Then combine overlapping or "book-ended" features
bedtools slop -i ${GENE_ANN} -g ${REF_FAI} -b ${RESIZE_FRACTION} -pct | sort -k1,1 -k2,2n | bedtools merge -i - > ${OUT_DIR}/${ann_prefix}.extended_${RESIZE_FRACTION}.merged.bed

# Prepare output file prefix
low_complexity_prefix=$(basename ${LOW_COMPLEXITY} .bed)
# Create low complexity regions that exclude 10% up or downstream of gene annotations
bedtools subtract -a ${LOW_COMPLEXITY} -b ${OUT_DIR}/${ann_prefix}.extended_${RESIZE_FRACTION}.merged.bed > ${OUT_DIR}/${low_complexity_prefix}.subtracted_gene_ann.bed

# Concatenate and merge uncallable regions into a single file
cat ${REF_Ns_BED} ${REPEAT_ANN} ${HIGH_COPY_BED} | grep -vw "chrUn" | cut -f 1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > ${OUT_DIR}/${OUT_PREFIX}.nochrUn.bed
# Add in low complexity regions excluding overlap with gene annotations
cat ${REF_Ns_BED} ${REPEAT_ANN} ${HIGH_COPY_BED} ${OUT_DIR}/${low_complexity_prefix}.subtracted_gene_ann.bed | grep -vw "chrUn" | cut -f 1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > ${OUT_DIR}/${OUT_PREFIX}.low_complexity.nochrUn.bed
