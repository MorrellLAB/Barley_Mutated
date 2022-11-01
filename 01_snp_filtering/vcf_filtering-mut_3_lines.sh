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
module load gatk/4.1.2
# Required for vcf-annotate that is part of vcftools
export PERL5LIB=$PERL5LIB:/panfs/jay/groups/9/morrellp/public/Software/vcftools_ML-0.1.16/share/perl5

# User provided input arguments
#RAW_VCF=/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated/mut_3_lines_sorted.vcf
# M01
DEL_VCF_M01="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M01-3-3/outs/M01-3-3_dels.vcf.gz"
LSV_VCF_M01="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/fix_ref_mismatch/M01-3-3_large_svs_noRefMismatch.vcf.gz"
PHV_VCF_M01="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/complex_variants/M01-3-3_phased_variants_noComplex.vcf.gz"
# M20
DEL_VCF_M20="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M20-2-2/outs/M20-2-2_dels.vcf.gz"
LSV_VCF_M20="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/fix_ref_mismatch/M20-2-2_large_svs_noRefMismatch.vcf.gz"
PHV_VCF_M20="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/complex_variants/M20-2-2_phased_variants_noComplex.vcf.gz"
# M29
DEL_VCF_M29="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M29-2-2/outs/M29-2-2_dels.vcf.gz"
LSV_VCF_M29="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/fix_ref_mismatch/M29-2-2_large_svs_noRefMismatch.vcf.gz"
PHV_VCF_M29="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/complex_variants/M29-2-2_phased_variants_noComplex.vcf.gz"
# Full filepath to output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered"
# Output file prefix
PREFIX_M01="M01-3-3"
PREFIX_M20="M20-2-2"
PREFIX_M29="M29-2-2"
PREFIX="mut_3_lines"

# Allele balance filter, minimum and maximum cutoff
# AB filter applied to phased variants only because on this VCF has the AD format field
MIN_AB="0.30"
MAX_AB="0.70"

VCFTOOLS_CUSTOM_FILTER="~/GitHub/Barley_Mutated/01_snp_filtering/filters.txt"

# List of regions where REF has stretches of N's
REF_Ns_BED="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/stretches_of_Ns/Barley_MorexV3_pseudomolecules_parts_missing.bed"
# Repeat annotations
REPEAT_ANN="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/annotation/HvulgareMorex_702_V3.repeatmasked_assembly_V3.parts.gff3"
# High copy regions (e.g., chloroplasts, mitochondria, rDNA repeats, centromere repeats, etc.)
HIGH_COPY_BED="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/high_copy_regions/Morex_v3_high_copy_uniq.parts.bed"

# BED file containing sites that differ between 10x morex-sample2 and Morex reference
# SNPs and Indels called in the phased variants VCF
PHV_MOREX_DIFFS_SNPs="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/morex-sample2_phased_variants-snps.DPfilt.noRepeatOverlap.noRefNs.diffs_from_ref.bed"
PHV_MOREX_DIFFS_INDELS="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/morex-sample2_phased_variants-indels.DPfilt.noRepeatOverlap.noRefNs.diffs_from_ref.bed"

# Known SNPs
SNPs_BOPA="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k/bopa_idt95_noRescuedSNPs_partsRef.vcf"
SNPs_9k="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k/9k_idt95_noRescuedSNPs_partsRef.vcf"
SNPs_50k="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k/50k_idt95_noRescuedSNPs_partsRef.vcf"

#----------------------------------------
# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR} ${OUT_DIR}/Intermediates

# Pass 1 filtering
# Filter out sites using 10x Genomics custom filters and exclude chrUn
function pass1_filtering() {
    local del_vcf="$1"
    local large_svs_vcf="$2"
    local phased_var_vcf="$3"
    local out_dir="$4"
    local sample_prefix="$5"
    ### DEL
    bcftools filter --targets "^chrUn" -e 'FILTER=="LOWQ"' ${del_vcf} -O z -o ${out_dir}/${sample_prefix}_dels.10xCustomFilt.vcf.gz
    # Separate BNDs from other variant types
    # File of other SVs, excluded BNDs
    bcftools view -e 'INFO/SVTYPE=="BND" | INFO/SVTYPE2=="BND"' ${out_dir}/${sample_prefix}_dels.10xCustomFilt.vcf.gz -O z -o ${out_dir}/${sample_prefix}_dels.10xCustomFilt.noBND.vcf.gz
    # File of only BND variants
    bcftools view -i 'INFO/SVTYPE=="BND" | INFO/SVTYPE2=="BND"' ${out_dir}/${sample_prefix}_dels.10xCustomFilt.vcf.gz -O z -o ${out_dir}/${sample_prefix}_dels.10xCustomFilt.BND_only.vcf.gz
    ### Large SVs
    bcftools filter --targets "^chrUn" -e 'FILTER=="LOWQ"' ${large_svs_vcf} | bcftools filter -e 'ALT~"chrUn"' -O z -o ${out_dir}/${sample_prefix}_large_svs.10xCustomFilt.vcf.gz
    # Separate BNDs from other variant types
    # File of other SVs, excluded BNDs
    bcftools view -e 'INFO/SVTYPE=="BND" | INFO/SVTYPE2=="BND" | INFO/SVTYPE=="UNK"' ${out_dir}/${sample_prefix}_large_svs.10xCustomFilt.vcf.gz -O z -o ${out_dir}/${sample_prefix}_large_svs.10xCustomFilt.noBND.vcf.gz
    # File of only BND variants
    bcftools view -i 'INFO/SVTYPE=="BND" | INFO/SVTYPE2=="BND" | INFO/SVTYPE=="UNK"' ${out_dir}/${sample_prefix}_large_svs.10xCustomFilt.vcf.gz -O z -o ${out_dir}/${sample_prefix}_large_svs.10xCustomFilt.BND_only.vcf.gz
    # Phased variants
    # Use 10x Genomics tuned filters and filter out sites with DP < 5
    bcftools filter --targets "^chrUn" -e 'FILTER=="10X_QUAL_FILTER" || FILTER=="10X_ALLELE_FRACTION_FILTER" || FILTER=="10X_PHASING_INCONSISTENT" || FILTER=="10X_HOMOPOLYMER_UNPHASED_INSERTION" || FILTER=="10X_RESCUED_MOLECULE_HIGH_DIVERSITY"' ${phased_var_vcf} | bcftools filter -e 'INFO/DP < 5' -O z -o ${out_dir}/${sample_prefix}_phased_variants.10xCustomFilt.vcf.gz
}

export -f pass1_filtering

function pass1_AB_filtering() {
    local phased_var_vcf="$1"
    local out_dir="$2"
    local sample_prefix="$3"
    local min_ab="$4"
    local max_ab="$5"
    # For phased variants only (has AD format field) that are heterozygotes, filter by allelic balance where cutoffs are set above
    # AB defined the same way as Pedersen et al. 2021: alt/(ref+alt)
    # Use bcftools to set genotypes to missing based on cutoffs
    #   -i in +setGT means if GT ann meet condition, set to missing
    #   FMT/AD[0:1] means first sample, second AD value
    # Then remove sites where single sample has been set to missing
    # This way allows us to build up the command and check if our expressions work as expected
    bcftools +setGT ${phased_var_vcf} -- -t q -n "." -i "(GT='het' & (FMT/AD[0:1])/(FMT/AD[0:0]+FMT/AD[0:1])<${min_ab}) | (GT='het' & (FMT/AD[0:1])/(FMT/AD[0:0]+FMT/AD[0:1])>${max_ab})" | bcftools filter -e 'GT="mis"' -O z -o ${out_dir}/${sample_prefix}_phased_variants.10xCustomFilt.ABfilt.vcf.gz
    # Index vcf
    tabix -p vcf ${out_dir}/${sample_prefix}_phased_variants.10xCustomFilt.ABfilt.vcf.gz
}

export -f pass1_AB_filtering

function pass2_filtering() {
    local del_vcf="$1"
    local large_svs_vcf="$2"
    local phased_var_vcf="$3"
    local out_dir="$4"
    local sample_prefix="$5"
    local repeat_ann="$6"
    local high_copy_bed="$7"
    local ref_ns_bed="$8"
    # Remove SVs that overlap with repeat annotated regions, high copy regions, and that overlap with stretches of Ns
    # DEL
    bedtools intersect -wa -v -header -a ${del_vcf} -b ${repeat_ann} ${high_copy_bed} ${ref_ns_bed} | bgzip > ${out_dir}/${sample_prefix}_dels.noRepeatOverlap.noRefNs.vcf.gz
    # Large SVs
    bedtools intersect -wa -v -header -a ${large_svs_vcf} -b ${repeat_ann} ${high_copy_bed} ${ref_ns_bed} | bgzip > ${out_dir}/${sample_prefix}_large_svs.noRepeatOverlap.noRefNs.vcf.gz
    # Phased variants
    bedtools intersect -wa -v -header -a ${phased_var_vcf} -b ${repeat_ann} ${high_copy_bed} ${ref_ns_bed} | bgzip > ${out_dir}/${sample_prefix}_phased_variants.noRepeatOverlap.noRefNs.vcf.gz
    # Index vcfs
    tabix -p vcf ${out_dir}/${sample_prefix}_dels.noRepeatOverlap.noRefNs.vcf.gz
    tabix -p vcf ${out_dir}/${sample_prefix}_large_svs.noRepeatOverlap.noRefNs.vcf.gz
    tabix -p vcf ${out_dir}/${sample_prefix}_phased_variants.noRepeatOverlap.noRefNs.vcf.gz
}

export -f pass2_filtering

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
        -O "${out_dir}/${out_bn}.INDELs.vcf.gz"
}

export -f separate_snps_and_indels

# Pass1 filtering: 10x Genomics custom filters and separating BND variants
pass1_filtering ${DEL_VCF_M01} ${LSV_VCF_M01} ${PHV_VCF_M01} ${OUT_DIR} ${PREFIX_M01}
pass1_filtering ${DEL_VCF_M20} ${LSV_VCF_M20} ${PHV_VCF_M20} ${OUT_DIR} ${PREFIX_M20}
pass1_filtering ${DEL_VCF_M29} ${LSV_VCF_M29} ${PHV_VCF_M29} ${OUT_DIR} ${PREFIX_M29}

# For phased variants only, filter on AB
# This takes the ${OUT_DIR}/${PREFIX}_phased_variants.10xCustomFilt.vcf.gz as input for each respective sample
pass1_AB_filtering ${OUT_DIR}/${PREFIX_M01}_phased_variants.10xCustomFilt.vcf.gz ${OUT_DIR} ${PREFIX_M01} ${MIN_AB} ${MAX_AB}
pass1_AB_filtering ${OUT_DIR}/${PREFIX_M20}_phased_variants.10xCustomFilt.vcf.gz ${OUT_DIR} ${PREFIX_M20} ${MIN_AB} ${MAX_AB}
pass1_AB_filtering ${OUT_DIR}/${PREFIX_M29}_phased_variants.10xCustomFilt.vcf.gz ${OUT_DIR} ${PREFIX_M29} ${MIN_AB} ${MAX_AB}

# Pass2 filtering
# Remove SVs that overlap with repeat annotated regions and that overlap with stretches of Ns
# M01
pass2_filtering \
    ${OUT_DIR}/${PREFIX_M01}_dels.10xCustomFilt.noBND.vcf.gz \
    ${OUT_DIR}/${PREFIX_M01}_large_svs.10xCustomFilt.noBND.vcf.gz \
    ${OUT_DIR}/${PREFIX_M01}_phased_variants.10xCustomFilt.ABfilt.vcf.gz \
    ${OUT_DIR} \
    ${PREFIX_M01} \
    ${REPEAT_ANN} \
    ${HIGH_COPY_BED} \
    ${REF_Ns_BED}
# M20
pass2_filtering \
    ${OUT_DIR}/${PREFIX_M20}_dels.10xCustomFilt.noBND.vcf.gz \
    ${OUT_DIR}/${PREFIX_M20}_large_svs.10xCustomFilt.noBND.vcf.gz \
    ${OUT_DIR}/${PREFIX_M20}_phased_variants.10xCustomFilt.ABfilt.vcf.gz \
    ${OUT_DIR} \
    ${PREFIX_M20} \
    ${REPEAT_ANN} \
    ${HIGH_COPY_BED} \
    ${REF_Ns_BED}
# M29
pass2_filtering \
    ${OUT_DIR}/${PREFIX_M29}_dels.10xCustomFilt.noBND.vcf.gz \
    ${OUT_DIR}/${PREFIX_M29}_large_svs.10xCustomFilt.noBND.vcf.gz \
    ${OUT_DIR}/${PREFIX_M29}_phased_variants.10xCustomFilt.ABfilt.vcf.gz \
    ${OUT_DIR} \
    ${PREFIX_M29} \
    ${REPEAT_ANN} \
    ${HIGH_COPY_BED} \
    ${REF_Ns_BED}

# Merge and concat files at this stage to explore visually in IGV
# DELs
bcftools merge \
    ${OUT_DIR}/${PREFIX_M01}_dels.noRepeatOverlap.noRefNs.vcf.gz \
    ${OUT_DIR}/${PREFIX_M20}_dels.noRepeatOverlap.noRefNs.vcf.gz \
    ${OUT_DIR}/${PREFIX_M29}_dels.noRepeatOverlap.noRefNs.vcf.gz \
    -O z -o ${OUT_DIR}/${PREFIX}_dels_merged.noRepeatOverlap.noRefNs.vcf.gz
# Large SVs
bcftools merge \
    ${OUT_DIR}/${PREFIX_M01}_large_svs.noRepeatOverlap.noRefNs.vcf.gz \
    ${OUT_DIR}/${PREFIX_M20}_large_svs.noRepeatOverlap.noRefNs.vcf.gz \
    ${OUT_DIR}/${PREFIX_M29}_large_svs.noRepeatOverlap.noRefNs.vcf.gz \
    -O z -o ${OUT_DIR}/${PREFIX}_large_svs_merged.noRepeatOverlap.noRefNs.vcf.gz
# Phased variants
bcftools merge \
    ${OUT_DIR}/${PREFIX_M01}_phased_variants.noRepeatOverlap.noRefNs.vcf.gz \
    ${OUT_DIR}/${PREFIX_M20}_phased_variants.noRepeatOverlap.noRefNs.vcf.gz \
    ${OUT_DIR}/${PREFIX_M29}_phased_variants.noRepeatOverlap.noRefNs.vcf.gz \
    -O z -o ${OUT_DIR}/${PREFIX}_phased_variants_merged.noRepeatOverlap.noRefNs.vcf.gz
# Index vcfs
tabix -p vcf ${OUT_DIR}/${PREFIX}_phased_variants_merged.noRepeatOverlap.noRefNs.vcf.gz
tabix -p vcf ${OUT_DIR}/${PREFIX}_dels_merged.noRepeatOverlap.noRefNs.vcf.gz
tabix -p vcf ${OUT_DIR}/${PREFIX}_large_svs_merged.noRepeatOverlap.noRefNs.vcf.gz
# Concatenate files for IGV exploration
bcftools concat --allow-overlaps \
    ${OUT_DIR}/${PREFIX}_phased_variants_merged.noRepeatOverlap.noRefNs.vcf.gz \
    ${OUT_DIR}/${PREFIX}_dels_merged.noRepeatOverlap.noRefNs.vcf.gz \
    ${OUT_DIR}/${PREFIX}_large_svs_merged.noRepeatOverlap.noRefNs.vcf.gz \
    | bcftools sort -O v -o ${OUT_DIR}/${PREFIX}_all_var_filt_concat.vcf

# Separate SNPs and indels in phased variants VCF
separate_snps_and_indels ${OUT_DIR}/${PREFIX_M01}_phased_variants.noRepeatOverlap.noRefNs.vcf.gz ${OUT_DIR}
separate_snps_and_indels ${OUT_DIR}/${PREFIX_M20}_phased_variants.noRepeatOverlap.noRefNs.vcf.gz ${OUT_DIR}
separate_snps_and_indels ${OUT_DIR}/${PREFIX_M29}_phased_variants.noRepeatOverlap.noRefNs.vcf.gz ${OUT_DIR}

# Exclude sites that differ between 10x Morex and Morex reference
# Also exclude sites that are known variants
# M01
bedtools intersect -v -header -a ${OUT_DIR}/${PREFIX_M01}_phased_variants.noRepeatOverlap.noRefNs.SNPs.vcf.gz -b ${PHV_MOREX_DIFFS_SNPs} ${SNPs_BOPA} ${SNPs_9k} ${SNPs_50k} | bgzip > ${OUT_DIR}/${PREFIX_M01}_phased_variants.noRepeatOverlap.noRefNs.SNPs.noMorexDiffs.vcf.gz
tabix -p vcf ${OUT_DIR}/${PREFIX_M01}_phased_variants.noRepeatOverlap.noRefNs.SNPs.noMorexDiffs.vcf.gz

bedtools intersect -v -header -a ${OUT_DIR}/${PREFIX_M01}_phased_variants.noRepeatOverlap.noRefNs.INDELs.vcf.gz -b ${PHV_MOREX_DIFFS_INDELS} | bgzip > ${OUT_DIR}/${PREFIX_M01}_phased_variants.noRepeatOverlap.noRefNs.INDELs.noMorexDiffs.vcf.gz
tabix -p vcf ${OUT_DIR}/${PREFIX_M01}_phased_variants.noRepeatOverlap.noRefNs.INDELs.noMorexDiffs.vcf.gz

# M20
bedtools intersect -v -header -a ${OUT_DIR}/${PREFIX_M20}_phased_variants.noRepeatOverlap.noRefNs.SNPs.vcf.gz -b ${PHV_MOREX_DIFFS_SNPs} ${SNPs_BOPA} ${SNPs_9k} ${SNPs_50k} | bgzip > ${OUT_DIR}/${PREFIX_M20}_phased_variants.noRepeatOverlap.noRefNs.SNPs.noMorexDiffs.vcf.gz
tabix -p vcf ${OUT_DIR}/${PREFIX_M20}_phased_variants.noRepeatOverlap.noRefNs.SNPs.noMorexDiffs.vcf.gz

bedtools intersect -v -header -a ${OUT_DIR}/${PREFIX_M20}_phased_variants.noRepeatOverlap.noRefNs.INDELs.vcf.gz -b ${PHV_MOREX_DIFFS_INDELS} | bgzip > ${OUT_DIR}/${PREFIX_M20}_phased_variants.noRepeatOverlap.noRefNs.INDELs.noMorexDiffs.vcf.gz
tabix -p vcf ${OUT_DIR}/${PREFIX_M20}_phased_variants.noRepeatOverlap.noRefNs.INDELs.noMorexDiffs.vcf.gz
# M29
bedtools intersect -v -header -a ${OUT_DIR}/${PREFIX_M29}_phased_variants.noRepeatOverlap.noRefNs.SNPs.vcf.gz -b ${PHV_MOREX_DIFFS_SNPs} ${SNPs_BOPA} ${SNPs_9k} ${SNPs_50k} | bgzip > ${OUT_DIR}/${PREFIX_M29}_phased_variants.noRepeatOverlap.noRefNs.SNPs.noMorexDiffs.vcf.gz
tabix -p vcf ${OUT_DIR}/${PREFIX_M29}_phased_variants.noRepeatOverlap.noRefNs.SNPs.noMorexDiffs.vcf.gz

bedtools intersect -v -header -a ${OUT_DIR}/${PREFIX_M29}_phased_variants.noRepeatOverlap.noRefNs.INDELs.vcf.gz -b ${PHV_MOREX_DIFFS_INDELS} | bgzip > ${OUT_DIR}/${PREFIX_M29}_phased_variants.noRepeatOverlap.noRefNs.INDELs.noMorexDiffs.vcf.gz
tabix -p vcf ${OUT_DIR}/${PREFIX_M29}_phased_variants.noRepeatOverlap.noRefNs.INDELs.noMorexDiffs.vcf.gz

# Cleanup, move intermediate filtering steps output to sub directory
# phased variants VCFs
mv ${OUT_DIR}/${PREFIX_M01}_phased_variants.10xCustomFilt.vcf.gz* \
    ${OUT_DIR}/${PREFIX_M20}_phased_variants.10xCustomFilt.vcf.gz* \
    ${OUT_DIR}/${PREFIX_M29}_phased_variants.10xCustomFilt.vcf.gz* \
    ${OUT_DIR}/${PREFIX_M01}_phased_variants.10xCustomFilt.ABfilt.vcf.gz* \
    ${OUT_DIR}/${PREFIX_M20}_phased_variants.10xCustomFilt.ABfilt.vcf.gz* \
    ${OUT_DIR}/${PREFIX_M29}_phased_variants.10xCustomFilt.ABfilt.vcf.gz* \
    ${OUT_DIR}/Intermediates
