#!/bin/bash

set -e
set -o pipefail

# Filter 10x Genomics VCF files for dels (larger deletions) and larger SVs and identify differences between 10x Morex and Morex reference
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
# VCF files: dels.vcf.gz, large_svs.vcf.gz
DEL_VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/outs/morex-sample2_dels.vcf.gz"
# Large SVs VCF: no chrUn, no DUP or UNK variant types, and no reference mismatches
LSV_VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/outs/morex-sample2_large_svs.calls.nochrUn.noDUP-UNK.vcf"
# Where do we output our files?
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered"
# Output file prefix
PREFIX="morex-sample2"

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

# Prepare large SVs VCF before filtering
LSV_DIR=$(dirname ${LSV_VCF})
if [[ "${LSV_VCF}" == *"vcf.gz" ]]; then
    LSV_BN=$(basename ${LSV_VCF} .vcf.gz)
elif [[ "${LSV_VCF}" == *"vcf" ]]; then
    LSV_BN=$(basename ${LSV_VCF} .vcf)
fi
# Separate: INV, DEL, and BND
# These will have different filtering criteria
bcftools view -i 'INFO/SVTYPE="DEL"' ${LSV_VCF} -O v -o ${LSV_DIR}/${LSV_BN}.DEL.vcf
bcftools view -i 'INFO/SVTYPE="INV"' ${LSV_VCF} -O v -o ${LSV_DIR}/${LSV_BN}.INV.vcf
bcftools view -i 'INFO/SVTYPE="BND" | INFO/SVTYPE2="BND"' ${LSV_VCF} -O v -o ${LSV_DIR}/${LSV_BN}.BND.vcf

# First pass filtering
# Filter out sites using 10x Genomics custom filters
# See here for info on filters:
#   https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/vcf
# Also remove sites where there are unplaced scaffolds (chrUn)
# DEL
bcftools filter --targets "^chrUn" -e 'FILTER=="LOWQ"' ${DEL_VCF} -O z -o ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.vcf.gz

# Large SVs
# Also exlude site if BND falls on chrUn (ALT column has chrUn), we can't say much about these
bcftools filter --targets "^chrUn" -e 'FILTER=="LOWQ"' ${LSV_DIR}/${LSV_BN}.DEL.vcf | bcftools filter -e 'ALT~"chrUn" | GT[*]="RR"' -O z -o ${OUT_DIR}/${PREFIX}_large_svs.DEL.10xCustomFilt.vcf.gz
bcftools filter --targets "^chrUn" -e 'FILTER=="LOWQ"' ${LSV_DIR}/${LSV_BN}.INV.vcf | bcftools filter -e 'ALT~"chrUn" | GT[*]="RR"' -O z -o ${OUT_DIR}/${PREFIX}_large_svs.INV.10xCustomFilt.vcf.gz
# At this point, look at large SVs INV in IGV and Loupe software
#   and decide if any of them look real. Here, we have less than a handful remaining after filtering
#   so it is manageable to look through in IGV and Loupe.

# Second pass filtering
# Remove sites that are homozygous reference, we want to identify sites that are differences from reference
# Also separate BND type SVs into separate files (BND only in dels.vcf.gz and large svs vcf)
# DEL
# Pull out variants without DP annotation that are NOT BND (SV breakend sets)
bcftools view -e 'INFO/SVTYPE=="BND" | INFO/SVTYPE2=="BND"' ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.vcf.gz | bcftools view -e 'GT[*]="RR"' -O z -o ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.noBND.vcf.gz
# Pull out BND (SV breakend sets)
#   These need to be dealt with differently when converting to BED format
bcftools view -i 'INFO/SVTYPE=="BND" | INFO/SVTYPE2=="BND"' ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.vcf.gz | bcftools view -e 'GT[*]="RR"' -O z -o ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.BND_only.vcf.gz
# Index vcf
bcftools index ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.noBND.vcf.gz
bcftools index ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.BND_only.vcf.gz

# Third pass filtering
# Remove SVs that overlap with uncallable regions relevant to SNPs and 1bp indels
# DELs (non BND types)
bedtools intersect -wa -v -header -a ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.noBND.vcf.gz -b ${REPEAT_ANN} ${HIGH_COPY_BED} ${REF_Ns_BED} ${LOW_COMPLEXITY} | bgzip > ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.noBND.callable.vcf.gz
# There are only ~208 larger DELs, so we'll visualize these using Samplot/SV-plaudit
#   and come up with a "pass" DELs list

# Large SVs (non BND types)
# Note: 0 variants remaining after filtering, so no need to include in downstream steps
bedtools intersect -wa -v -header -a ${OUT_DIR}/${PREFIX}_large_svs.DEL.10xCustomFilt.vcf.gz -b ${REPEAT_ANN} ${HIGH_COPY_BED} ${REF_Ns_BED} | bgzip > ${OUT_DIR}/${PREFIX}_large_svs.DEL.10xCustomFilt.noBND.callable.vcf.gz
# Index vcfs
tabix -p vcf ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.noBND.callable.vcf.gz
tabix -p vcf ${OUT_DIR}/${PREFIX}_large_svs.DEL.10xCustomFilt.noBND.callable.vcf.gz

# Convert diffs from ref VCF to BED
# Use a custom script for 10x Genomics specific dels, dups, and SVs formatting
# These files will be used as exclusion lists
# Custom script for 10x Genomics specific format that are NOT BND (breakend sets)
# Only for dels.vcf.gz set since large SVs had 0 SVs remain after filtering
# Remember, these are the SV-plaudit curated DELs
vcf_10x_genomics_to_bed.py ${SCORED_VCF} > ${OUT_DIR}/diffs_from_ref/${PREFIX}_dels_diffs_from_ref.noBND.bed

# Merge BED files
# SNPs that fall in morex larger SV regions will be counted as part of the bigger diffs from ref region
#cat ${OUT_DIR}/diffs_from_ref/${PREFIX}_dels_diffs_from_ref.noBND.bed ${OUT_DIR}/diffs_from_ref/${PREFIX}_phased_variants-snps.DPfilt.callable.diffs_from_ref.bed ${OUT_DIR}/diffs_from_ref/${PREFIX}_phased_variants-indels.DPfilt.callable.diffs_from_ref.bed | sort -k1,1 -k2,2n | bedtools merge -i - > ${OUT_DIR}/diffs_from_ref/${PREFIX}_combined.diffs_from_ref.bed

# Cleanup, move intermediate filtering steps output to sub directory
# larger indels vcfs
mv ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.vcf.gz* \
    ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.noBND.vcf.gz* \
    ${OUT_DIR}/Intermediates
# larger SVs vcfs
mv ${OUT_DIR}/${PREFIX}_large_svs.DEL.10xCustomFilt.vcf.gz* \
    ${OUT_DIR}/Intermediates
