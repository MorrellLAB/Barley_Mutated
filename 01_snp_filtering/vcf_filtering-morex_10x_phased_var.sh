#!/bin/bash

set -e
set -o pipefail

# Filter 10x Genomics VCF file (includes SNPs, dels, and SVs) and identify differences between 10x Morex and Morex reference
#   This script outputs a VCF with differences from reference

# Dependencies
module load bcftools/1.9
module load bedtools/2.29.2
module load python3/3.8.3_anaconda2020.07_mamba
module load vcflib_ML/1.0.0_rc2
module load htslib/1.9
module load gatk/4.1.2
# Export path to directory containing custom script for converting 10x Genomics VCF
#   containing dels, dups, and SVs where the end position is in the INFO field
export PATH=${PATH}:/panfs/jay/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/01_snp_filtering

# User provided input arguments
# VCF files: dels.vcf.gz, large_svs.vcf.gz, phased_variants.vcf.gz
# Phased variants VCF
PHV_VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/outs/morex-sample2_phased_variants.vcf.gz"
# Where do we output our files?
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered"
# Output file prefix
PREFIX="morex-sample2"

# Allele balance filter, minimum and maximum cutoff
# AB filter applied to phased variants only because on this VCF has the AD format field
MIN_AB="0.30"
MAX_AB="0.70"

# Minimum and maximum per sample DP
MIN_DP="5"
MAX_DP="78"

# List of regions where REF has stretches of N's
REF_Ns_BED="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/stretches_of_Ns/Barley_MorexV3_pseudomolecules_parts_missing.bed"
# Repeat annotations
REPEAT_ANN="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/annotation/HvulgareMorex_702_V3.repeatmasked_assembly_V3.parts.gff3"
# High copy regions (e.g., chloroplasts, mitochondria, rDNA repeats, centromere repeats, etc.)
HIGH_COPY_BED="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/high_copy_regions/Morex_v3_high_copy_uniq.parts.bed"
# Low complexity regions generated from JGI's BBMask
LOW_COMPLEXITY="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/entropy_masked/Barley_MorexV3_pseudomolecules_parts.entropy_0.7_masked.bed"
# High diversity, >2% diversity in a 400bp window for morex-sample2
HIGH_DIV_BED="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/uncallable_regions/pixy_pi_400bp_win.gt0.02.bed"

# SV-Plaudit scored VCF (supports only)
# Note: This is generated from script in the subdirectory Samplot-Morex
SCORED_VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/samplot-morex_10x/morex-sample2_dels.10xCustomFilt.noBND.callable.supports.vcf"

#---------------------
# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR} ${OUT_DIR}/Intermediates ${OUT_DIR}/diffs_from_ref

# First pass filtering
# Filter out sites using 10x Genomics custom filters
# See here for info on filters:
#   https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/vcf
# Also remove sites where there are unplaced scaffolds (chrUn)
# Phased variants
bcftools filter --targets "^chrUn" -e 'FILTER=="10X_QUAL_FILTER" || FILTER=="10X_ALLELE_FRACTION_FILTER" || FILTER=="10X_PHASING_INCONSISTENT" || FILTER=="10X_HOMOPOLYMER_UNPHASED_INSERTION" || FILTER=="10X_RESCUED_MOLECULE_HIGH_DIVERSITY"' ${PHV_VCF} -O z -o ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.10xCustomFilt.vcf.gz
# For phased variants only (has AD format field) that are heterozygotes, filter by allelic balance where cutoffs are set above
# AB defined the same way as Pedersen et al. 2021: alt/(ref+alt)
# Use bcftools to set genotypes to missing based on cutoffs
#   -i in +setGT means if GT ann meet condition, set to missing
#   FMT/AD[0:1] means first sample, second AD value
# Then remove sites where single sample has been set to missing
# This way allows us to build up the command and check if our expressions work as expected
bcftools +setGT ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.10xCustomFilt.vcf.gz -- -t q -n "." -i "(GT='het' & (FMT/AD[0:1])/(FMT/AD[0:0]+FMT/AD[0:1])<${MIN_AB}) | (GT='het' & (FMT/AD[0:1])/(FMT/AD[0:0]+FMT/AD[0:1])>${MAX_AB})" | bcftools filter -e 'GT="mis"' -O z -o ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.10xCustomFilt.ABfilt.vcf.gz
# Index vcf
tabix -p vcf ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.10xCustomFilt.ABfilt.vcf.gz

# A) Separate variants with the DP annotation from ones that don't have the DP annotation
# This is unique to the 10x Genomics VCF format where phased variants (i.e., SNPs, some indels)
#   have a DP annotation, but the dels and SVs VCFs don't have the DP annotation
#   so we have to filter them separately.
# B) For variants with DP annotation, keep sites where DP > 5 and DP < 78
# In this case there is only one sample so the INFO/DP and FORMAT DP values when they both exist
#   are the same.
#   DP of 78 is approximately 2 sd away from the mean here
# C) Remove sites that are homozygous reference, we want to identify sites that are differences from reference
bcftools filter -i "INFO/DP>${MIN_DP} && INFO/DP<${MAX_DP}" ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.10xCustomFilt.ABfilt.vcf.gz | bcftools view -e 'GT[*]="RR"' -O z -o ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.DPfilt.vcf.gz
# Index vcf
tabix -p vcf ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.DPfilt.vcf.gz

# Remove SVs that overlap with uncallable regions relevant to SNPs and indels
# Phased variants
bedtools intersect -wa -v -header -a ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.DPfilt.vcf.gz -b ${REPEAT_ANN} ${HIGH_COPY_BED} ${REF_Ns_BED} ${HIGH_DIV_BED} | bgzip > ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.DPfilt.callable.vcf.gz
tabix -p vcf ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.DPfilt.callable.vcf.gz

# Separate biallelic and multiallelic variants
bcftools view -m2 -M2 ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.DPfilt.callable.vcf.gz -O z -o ${OUT_DIR}/${PREFIX}_phased_variants.DPfilt.callable.biallelic.vcf.gz
bcftools view -m3 ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.DPfilt.callable.vcf.gz -O z -o ${OUT_DIR}/${PREFIX} _phased_variants.DPfilt.callable.multiallelic.vcf.gz
# Index vcf
tabix -p vcf ${OUT_DIR}/${PREFIX}_phased_variants.DPfilt.callable.biallelic.vcf.gz

# For phased variants, separate SNPs and indels
# Select SNPs only from phased variants VCF
gatk SelectVariants \
    -V "${OUT_DIR}/${PREFIX}_phased_variants.DPfilt.callable.biallelic.vcf.gz" \
    -select-type SNP \
    -O "${OUT_DIR}/diffs_from_ref/${PREFIX}_phased_variants-snps.callable.biallelic.vcf.gz"
# Select indels only from phased variants VCF
gatk SelectVariants \
    -V "${OUT_DIR}/${PREFIX}_phased_variants.DPfilt.callable.biallelic.vcf.gz" \
    -select-type INDEL \
    -O "${OUT_DIR}/diffs_from_ref/${PREFIX}_phased_variants-indels.callable.biallelic.vcf.gz"

# Separate 1 bp indels from >1bp indels
bcftools view -i 'STRLEN(REF) <= 2 & STRLEN(ALT) <= 2' ${OUT_DIR}/diffs_from_ref/${PREFIX}_phased_variants-indels.callable.biallelic.vcf.gz -O z -o ${OUT_DIR}/diffs_from_ref/${PREFIX}_phased_variants-indels.callable.biallelic.1bp.vcf.gz
bcftools view -i 'STRLEN(REF) > 2 | STRLEN(ALT) > 2' ${OUT_DIR}/diffs_from_ref/${PREFIX}_phased_variants-indels.callable.biallelic.vcf.gz -O z -o ${OUT_DIR}/diffs_from_ref/${PREFIX}_phased_variants-indels.callable.biallelic.gt1bp.vcf.gz

# Convert diffs from ref VCF to BED using bedops tool for typical VCF format (phased variants)
# and use a custom script for 10x Genomics specific dels, dups, and SVs formatting
# These files will be used as exclusion lists
# Use vcflib to convert vcf to bed
# snps
zcat ${OUT_DIR}/diffs_from_ref/${PREFIX}_phased_variants-snps.callable.biallelic.vcf.gz | vcf2bed.py - | cut -f 1,2,3 > ${OUT_DIR}/diffs_from_ref/${PREFIX}_phased_variants-snps.callable.biallelic.diffs_from_ref.bed
# indels (1 bp from phased variants VCF)
#zcat ${OUT_DIR}/diffs_from_ref/${PREFIX}_phased_variants-indels.DPfilt.callable.vcf.gz | vcf2bed.py - | cut -f 1,2,3 > ${OUT_DIR}/diffs_from_ref/${PREFIX}_phased_variants-indels.DPfilt.callable.diffs_from_ref.bed
zcat ${OUT_DIR}/diffs_from_ref/${PREFIX}_phased_variants-indels.callable.biallelic.1bp.vcf.gz | vcf2bed.py - | cut -f 1,2,3 > ${OUT_DIR}/diffs_from_ref/${PREFIX}_phased_variants-indels.callable.biallelic.1bp.diffs_from_ref.bed
zcat ${OUT_DIR}/diffs_from_ref/${PREFIX}_phased_variants-indels.callable.biallelic.gt1bp.vcf.gz | vcf2bed.py - | cut -f 1,2,3 > ${OUT_DIR}/diffs_from_ref/${PREFIX}_phased_variants-indels.callable.biallelic.gt1bp.diffs_from_ref.bed
