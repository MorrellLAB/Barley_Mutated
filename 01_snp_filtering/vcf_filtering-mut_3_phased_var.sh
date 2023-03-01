#!/bin/bash

set -e
set -o pipefail

# Filter 10x Genomics VCF file containing 3 mutated lines (includes SNPs, dels, and SVs)
# Output VCF file containing singletons

# Dependencies
module load bcftools/1.9
module load htslib/1.9
module load bedtools/2.29.2
module load vcftools_ML/0.1.16
module load java/openjdk-8_202
module load gatk/4.1.2
# Required for vcf-annotate that is part of vcftools
export PERL5LIB=$PERL5LIB:/panfs/jay/groups/9/morrellp/public/Software/vcftools_ML-0.1.16/share/perl5

# User provided input arguments
# List of 10x Genomics 3 mutated lines phased variants VCF
PHV_VCF_LIST="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/phased_var_mut_only_vcf_list.txt"
# Full filepath to output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered"
# Output file prefix
PREFIX="mut_3_lines"

# Allele balance filter, minimum and maximum cutoff
# AB filter applied to phased variants only because on this VCF has the AD format field
MIN_AB="0.30"
MAX_AB="0.70"

# Min DP per sample threshold
MIN_DP="5"

# Uncallable BED file includes: REF has stretches of N's, repeat annotations, and high copy regions
UNCALLABLE_BED="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/uncallable_regions/morex_v3_combined_uncallable.nochrUn.bed"
# Uncallable BED file that includes the above plus low complexity regions
UNCALLABLE_LOW_COMPLEXITY_BED="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/uncallable_regions/morex_v3_combined_uncallable.low_complexity.nochrUn.bed"
# High diversity, >2% diversity in a 400bp window for morex-sample2
HIGH_DIV_BED="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/uncallable_regions/pixy_pi_400bp_win.gt0.02.bed"

# BED file containing sites that differ between 10x morex-sample2 and Morex reference
# SNPs and Indels called in the phased variants VCF
PHV_MOREX_DIFFS_SNPs="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/diffs_from_ref/morex-sample2_phased_variants-snps.callable.biallelic.diffs_from_ref.bed"
PHV_MOREX_DIFFS_1bp_INDELS="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/diffs_from_ref/morex-sample2_phased_variants-indels.callable.biallelic.1bp.diffs_from_ref.bed"
PHV_MOREX_DIFFS_gt1bp_INDELS="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/diffs_from_ref/morex-sample2_phased_variants-indels.callable.biallelic.gt1bp.diffs_from_ref.bed"

# Known SNPs
SNPs_BOPA="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k/bopa_idt95_noRescuedSNPs_partsRef.vcf"
SNPs_9k="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k/9k_idt95_noRescuedSNPs_partsRef.vcf"
SNPs_50k="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k/50k_idt95_noRescuedSNPs_partsRef.vcf"

#----------------------------------------
function separate_snps_and_indels() {
    local vcf="$1"
    local out_dir="$2"
    # For phased variants VCF, separate SNPs and indels
    # Prepare output file basename
    if [[ "${vcf}" == *"gz"* ]]; then
        # We are working with gzipped VCF file
        out_bn=$(basename ${vcf} .vcf.gz)
    else
        # We are working with uncompressed vcf
        out_bn=$(basename ${vcf} .vcf)
    fi
    # Select SNPs only from phased variants VCF
    gatk SelectVariants \
        -V "${vcf}" \
        -select-type SNP \
        -O "${out_dir}/${out_bn}.SNPs.vcf.gz"
    # Select indels only from phased variants VCF
    gatk SelectVariants \
        -V "${vcf}" \
        -select-type INDEL \
        -O "${out_dir}/${out_bn}.indels.vcf.gz"
}

export -f separate_snps_and_indels

# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR} ${OUT_DIR}/Intermediates

# Merge VCFs containing one sample each into one file
# and exclude chrUn
bcftools merge --file-list ${PHV_VCF_LIST} | bcftools filter --targets "^chrUn" -O z -o ${OUT_DIR}/${PREFIX}_phased_variants_merged.nochrUn.vcf.gz
tabix -p vcf ${OUT_DIR}/${PREFIX}_phased_variants_merged.nochrUn.vcf.gz

# Phased variants
# Use 10x Genomics tuned filters
# Filter to variants private to each mutated line (variants present in only a single sample)
bcftools filter -e 'FILTER=="10X_QUAL_FILTER" || FILTER=="10X_ALLELE_FRACTION_FILTER" || FILTER=="10X_PHASING_INCONSISTENT" || FILTER=="10X_HOMOPOLYMER_UNPHASED_INSERTION" || FILTER=="10X_RESCUED_MOLECULE_HIGH_DIVERSITY"' ${OUT_DIR}/${PREFIX}_phased_variants_merged.nochrUn.vcf.gz | bcftools view -i "COUNT(GT='alt')=1" -O z -o ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.10xCustomFilt.private.vcf.gz

#bcftools filter -e "INFO/DP < ${MIN_DP}" -O z -o ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.10xCustomFilt.vcf.gz

# bcftools +setGT ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.10xCustomFilt.biallelic.vcf.gz -- -t q -n "." -i "(GT='het' & (FMT/AD[*:1])/(FMT/AD[*:0]+FMT/AD[*:1])<${MIN_AB}) | (GT='het' & (FMT/AD[*:1])/(FMT/AD[*:0]+FMT/AD[*:1])>${MAX_AB})" | bcftools +setGT - -- -t q -n "." -i "FMT/DP<${MIN_DP}" | bcftools filter -e 'F_PASS(GT="mis") = 1.0' -O z -o ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.10xCustomFilt.biallelic.ABfilt.vcf.gz
# # Index vcf
# tabix -p vcf ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.10xCustomFilt.biallelic.ABfilt.vcf.gz

# Separate biallelic and multiallelic sites
# For phased variants only and biallelic variants (has AD format field) that are heterozygotes
# Set GT to missing for het that fail allelic balance threshold
# AB defined the same way as Pedersen et al. 2021: alt/(ref+alt)
# Use bcftools to set genotypes to missing based on cutoffs
#   -i in +setGT means if GT ann meet condition, set to missing
#   FMT/AD[0:1] means first sample, second AD value
# Also filter out sites with DP < MIN_DP
# biallelic
bcftools view -m2 -M2 ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.10xCustomFilt.private.vcf.gz | bcftools +setGT - -- -t q -n "." -i "(GT='het' & (FMT/AD[*:1])/(FMT/AD[*:0]+FMT/AD[*:1])<${MIN_AB}) | (GT='het' & (FMT/AD[*:1])/(FMT/AD[*:0]+FMT/AD[*:1])>${MAX_AB})" | bcftools +setGT - -- -t q -n "." -i "FMT/DP<${MIN_DP}" | bcftools filter -e 'F_PASS(GT="mis") = 1.0' -O z -o ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.10xCustomFilt.AB_DP.private.biallelic.vcf.gz
# multiallelic
bcftools view -m3 ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.10xCustomFilt.private.vcf.gz | bcftools +setGT - -- -t q -n "." -i "FMT/DP<${MIN_DP}" | bcftools filter -e 'F_PASS(GT="mis") = 1.0' -O z -o ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.10xCustomFilt.DP.private.multiallelic.vcf.gz

# Remove SVs that overlap with repeat annotated regions, high copy regions, and that overlap with stretches of Ns
# Phased variants
bedtools intersect -wa -v -header -a ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.10xCustomFilt.AB_DP.private.biallelic.vcf.gz -b ${UNCALLABLE_BED} ${HIGH_DIV_BED} | bgzip > ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.private.callable.vcf.gz
# Index vcf
tabix -p vcf ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.private.callable.vcf.gz

# Separate SNPs and indels in phased variants VCF

# Exclude sites that differ between 10x Morex and Morex reference
# Also exclude sites that are known variants
# SNPs and indels
bedtools intersect -v -header -a ${OUT_DIR}/Intermediates/${PREFIX}_phased_variants.private.callable.vcf.gz -b ${PHV_MOREX_DIFFS_SNPs} ${PHV_MOREX_DIFFS_1bp_INDELS} ${PHV_MOREX_DIFFS_gt1bp_INDELS} ${SNPs_BOPA} ${SNPs_9k} ${SNPs_50k} | bgzip > ${OUT_DIR}/${PREFIX}_phased_variants.private.callable.noMorexDiffs.vcf.gz
# Index vcf
tabix -p vcf ${OUT_DIR}/${PREFIX}_phased_variants.private.callable.noMorexDiffs.vcf.gz

# Separate SNPs and indels in phased variants VCF
separate_snps_and_indels ${OUT_DIR}/${PREFIX}_phased_variants.private.callable.noMorexDiffs.vcf.gz ${OUT_DIR}
