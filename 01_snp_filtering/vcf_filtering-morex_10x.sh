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
DEL_VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/outs/morex-sample2_dels.vcf.gz"
LSV_VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/fix_ref_mismatch/morex-sample2_large_svs_noRefMismatch.vcf.gz"
PHV_VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/complex_variants/morex-sample2_phased_variants_noComplex.vcf.gz"
# Where do we output our files?
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered"
# Output file prefix
PREFIX="morex-sample2"

# Allele balance filter, minimum and maximum cutoff
# AB filter applied to phased variants only because on this VCF has the AD format field
MIN_AB="0.30"
MAX_AB="0.70"

# List of regions where REF has stretches of N's
REF_Ns_BED="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/stretches_of_Ns/Barley_MorexV3_pseudomolecules_parts_missing.bed"
# Repeat annotations
REPEAT_ANN="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/annotation/HvulgareMorex_702_V3.repeatmasked_assembly_V3.parts.gff3"
# High copy regions (e.g., chloroplasts, mitochondria, rDNA repeats, centromere repeats, etc.)
HIGH_COPY_BED="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/high_copy_regions/Morex_v3_high_copy_uniq.parts.bed"

#---------------------
# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR} ${OUT_DIR}/Intermediates

# First pass filtering
# Filter out sites using 10x Genomics custom filters
# See here for info on filters:
#   https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/vcf
# Also remove sites where there are unplaced scaffolds (chrUn)
# DEL
bcftools filter --targets "^chrUn" -e 'FILTER=="LOWQ"' ${DEL_VCF} -O z -o ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.vcf.gz
# Large SVs
# Also exlude site if BND falls on chrUn (ALT column has chrUn), we can't say much about these
bcftools filter --targets "^chrUn" -e 'FILTER=="LOWQ"' ${LSV_VCF} | bcftools filter -e 'ALT~"chrUn"' -O z -o ${OUT_DIR}/${PREFIX}_large_svs.10xCustomFilt.vcf.gz
# Phased variants
bcftools filter --targets "^chrUn" -e 'FILTER=="10X_QUAL_FILTER" || FILTER=="10X_ALLELE_FRACTION_FILTER" || FILTER=="10X_PHASING_INCONSISTENT" || FILTER=="10X_HOMOPOLYMER_UNPHASED_INSERTION" || FILTER=="10X_RESCUED_MOLECULE_HIGH_DIVERSITY"' ${PHV_VCF} -O z -o ${OUT_DIR}/${PREFIX}_phased_variants.10xCustomFilt.vcf.gz
# For phased variants only (has AD format field) that are heterozygotes, filter by allelic balance where cutoffs are set above
# AB defined the same way as Pedersen et al. 2021: alt/(ref+alt)
# Use bcftools to set genotypes to missing based on cutoffs
#   -i in +setGT means if GT ann meet condition, set to missing
#   FMT/AD[0:1] means first sample, second AD value
# Then remove sites where single sample has been set to missing
# This way allows us to build up the command and check if our expressions work as expected
bcftools +setGT ${OUT_DIR}/${PREFIX}_phased_variants.10xCustomFilt.vcf.gz -- -t q -n "." -i "(GT='het' & (FMT/AD[0:1])/(FMT/AD[0:0]+FMT/AD[0:1])<${MIN_AB}) | (GT='het' & (FMT/AD[0:1])/(FMT/AD[0:0]+FMT/AD[0:1])>${MAX_AB})" | bcftools filter -e 'GT="mis"' -O z -o ${OUT_DIR}/${PREFIX}_phased_variants.10xCustomFilt.ABfilt.vcf.gz
# Index vcf
tabix -p vcf ${OUT_DIR}/${PREFIX}_phased_variants.10xCustomFilt.ABfilt.vcf.gz

# Second pass filtering
# A) Separate variants with the DP annotation from ones that don't have the DP annotation
# This is unique to the 10x Genomics VCF format where phased variants (i.e., SNPs, some indels)
#   have a DP annotation, but the dels and SVs VCFs don't have the DP annotation
#   so we have to filter them separately.
# B) For variants with DP annotation, keep sites where DP > 5 and DP < 78
# In this case there is only one sample so the INFO/DP and FORMAT DP values when they both exist
#   are the same.
#   DP of 78 is approximately 2 sd away from the mean here
# C) Remove sites that are homozygous reference, we want to identify sites that are differences from reference
# D) Also separate BND type SVs into separate files (BND only in dels.vcf.gz and large svs vcf)
# Phased variants
bcftools filter -i 'INFO/DP > 5 && INFO/DP < 78' ${OUT_DIR}/${PREFIX}_phased_variants.10xCustomFilt.ABfilt.vcf.gz | bcftools view -e 'GT[*]="RR"' -O z -o ${OUT_DIR}/${PREFIX}_phased_variants.DPfilt.vcf.gz
# Index vcf
bcftools index ${OUT_DIR}/${PREFIX}_phased_variants.DPfilt.vcf.gz

# DEL
# Pull out variants without DP annotation that are NOT BND (SV breakend sets)
bcftools view -e 'INFO/SVTYPE=="BND" | INFO/SVTYPE2=="BND"' ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.vcf.gz | bcftools view -e 'GT[*]="RR"' -O z -o ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.noBND.vcf.gz
# Pull out BND (SV breakend sets)
#   These need to be dealt with differently when converting to BED format
bcftools view -i 'INFO/SVTYPE=="BND" | INFO/SVTYPE2=="BND"' ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.vcf.gz | bcftools view -e 'GT[*]="RR"' -O z -o ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.BND_only.vcf.gz
# Index vcf
bcftools index ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.noBND.vcf.gz
bcftools index ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.BND_only.vcf.gz

# Large SVs
# Repeat these two steps for Large SVs
bcftools view -e 'INFO/SVTYPE=="BND" | INFO/SVTYPE2=="BND" | INFO/SVTYPE=="UNK"' ${OUT_DIR}/${PREFIX}_large_svs.10xCustomFilt.vcf.gz | bcftools view -e 'GT[*]="RR"' -O z -o ${OUT_DIR}/${PREFIX}_large_svs.10xCustomFilt.noBND.vcf.gz
bcftools view -i 'INFO/SVTYPE=="BND" | INFO/SVTYPE2=="BND" | INFO/SVTYPE=="UNK"' ${OUT_DIR}/${PREFIX}_large_svs.10xCustomFilt.vcf.gz | bcftools view -e 'GT[*]="RR"' -O z -o ${OUT_DIR}/${PREFIX}_large_svs.10xCustomFilt.BND_only.vcf.gz
# Index vcf
bcftools index ${OUT_DIR}/${PREFIX}_large_svs.10xCustomFilt.noBND.vcf.gz
bcftools index ${OUT_DIR}/${PREFIX}_large_svs.10xCustomFilt.BND_only.vcf.gz

# Third pass filtering
# Remove SVs that overlap with repeat annotated regions and that overlap with stretches of Ns
# Phased variants
bedtools intersect -wa -v -header -a ${OUT_DIR}/${PREFIX}_phased_variants.DPfilt.vcf.gz -b ${REPEAT_ANN} ${HIGH_COPY_BED} ${REF_Ns_BED} | bgzip > ${OUT_DIR}/${PREFIX}_phased_variants.DPfilt.noRepeatOverlap.noRefNs.vcf.gz
# DELs (non BND types)
bedtools intersect -wa -v -header -a ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.noBND.vcf.gz -b ${REPEAT_ANN} ${HIGH_COPY_BED} ${REF_Ns_BED} | bgzip > ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.noBND.noRepeatOverlap.noRefNs.vcf.gz
# Large SVs (non BND types)
bedtools intersect -wa -v -header -a ${OUT_DIR}/${PREFIX}_large_svs.10xCustomFilt.noBND.vcf.gz -b ${REPEAT_ANN} ${HIGH_COPY_BED} ${REF_Ns_BED} | bgzip > ${OUT_DIR}/${PREFIX}_large_svs.10xCustomFilt.noBND.noRepeatOverlap.noRefNs.vcf.gz
# Index vcfs
tabix -p vcf ${OUT_DIR}/${PREFIX}_phased_variants.DPfilt.noRepeatOverlap.noRefNs.vcf.gz
tabix -p vcf ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.noBND.noRepeatOverlap.noRefNs.vcf.gz
tabix -p vcf ${OUT_DIR}/${PREFIX}_large_svs.10xCustomFilt.noBND.noRepeatOverlap.noRefNs.vcf.gz

# Concatenate quality filtered VCFs for visualization with samplot, include BND variants
bcftools concat --allow-overlaps \
    ${OUT_DIR}/${PREFIX}_phased_variants.DPfilt.vcf.gz \
    ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.noBND.vcf.gz \
    ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.BND_only.vcf.gz \
    ${OUT_DIR}/${PREFIX}_large_svs.10xCustomFilt.noBND.vcf.gz \
    ${OUT_DIR}/${PREFIX}_large_svs.10xCustomFilt.BND_only.vcf.gz \
    | bcftools sort -O z -o ${OUT_DIR}/${PREFIX}_all_var_filt_concat.vcf.gz
# Index vcf
tabix -p vcf ${OUT_DIR}/${PREFIX}_all_var_filt_concat.vcf.gz
# Concatenate third pass filtering VCFs
bcftools concat --allow-overlaps \
    ${OUT_DIR}/${PREFIX}_phased_variants.DPfilt.noRepeatOverlap.noRefNs.vcf.gz \
    ${OUT_DIR}/${PREFIX}_dels.10xCustomFilt.noBND.noRepeatOverlap.noRefNs.vcf.gz \
    ${OUT_DIR}/${PREFIX}_large_svs.10xCustomFilt.noBND.noRepeatOverlap.noRefNs.vcf.gz \
    | bcftools sort -O z -o ${OUT_DIR}/${PREFIX}_all_var_filt_concat.noBND.noRepeatOverlap.noRefNs.vcf.gz
# Index vcf
tabix -p vcf ${OUT_DIR}/${PREFIX}_all_var_filt_concat.noBND.noRepeatOverlap.noRefNs.vcf.gz

# For phased variants, separate SNPs and indels
# Select SNPs only from phased variants VCF
gatk SelectVariants \
    -V "${OUT_DIR}/${PREFIX}_phased_variants.DPfilt.noRepeatOverlap.noRefNs.vcf.gz" \
    -select-type SNP \
    -O "${OUT_DIR}/${PREFIX}_phased_variants-snps.DPfilt.noRepeatOverlap.noRefNs.vcf.gz"
# Select indels only from phased variants VCF
gatk SelectVariants \
    -V "${OUT_DIR}/${PREFIX}_phased_variants.DPfilt.noRepeatOverlap.noRefNs.vcf.gz" \
    -select-type INDEL \
    -O "${OUT_DIR}/${PREFIX}_phased_variants-indels.DPfilt.noRepeatOverlap.noRefNs.vcf.gz"

# Convert diffs from ref VCF to BED using bedops tool for typical VCF format (phased variants)
# and use a custom script for 10x Genomics specific dels, dups, and SVs formatting
# These files will be used as exclusion lists
# Use vcflib to convert vcf to bed
# snps
zcat ${OUT_DIR}/${PREFIX}_phased_variants-snps.DPfilt.noRepeatOverlap.noRefNs.vcf.gz | vcf2bed.py - | cut -f 1,2,3 > ${OUT_DIR}/${PREFIX}_phased_variants-snps.DPfilt.noRepeatOverlap.noRefNs.diffs_from_ref.bed
# indels (1 bp from phased variants VCF)
zcat ${OUT_DIR}/${PREFIX}_phased_variants-indels.DPfilt.noRepeatOverlap.noRefNs.vcf.gz | vcf2bed.py - | cut -f 1,2,3 > ${OUT_DIR}/${PREFIX}_phased_variants-indels.DPfilt.noRepeatOverlap.noRefNs.diffs_from_ref.bed
# Custom script for 10x Genomics specific format that are NOT BND (breakend sets)
#vcf_10x_genomics_to_bed.py ${OUT_DIR}/${PREFIX}_filtered.noBND.vcf.gz > ${OUT_DIR}/${PREFIX}_diffs_from_ref.noBND.bed
# BND sites VCF to BEDPE
#bedtools intersect -wa -a ${BEDPE} -b ${OUT_DIR}/${PREFIX}_filtered.BND_only.vcf.gz | sort -uV -k1,3 > ${OUT_DIR}/${PREFIX}_diffs_from_ref.BND_only.bedpe

# Cleanup, move intermediate filtering steps output to sub directory
# phased variants VCFs
mv ${OUT_DIR}/${PREFIX}_phased_variants.10xCustomFilt.vcf.gz* \
    ${OUT_DIR}/${PREFIX}_phased_variants.10xCustomFilt.ABfilt.vcf.gz* \
    ${OUT_DIR}/${PREFIX}_phased_variants.DPfilt.vcf.gz* \
    ${OUT_DIR}/Intermediates
# larger indels vcfs

# larger SVs vcfs
