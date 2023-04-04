#!/bin/bash

set -e
set -o pipefail

# Exclude sites where there are differences between 10x Genomics morex-sample2 and Morex reference
# for the 8 WGS mutated lines

# Dependencies
module load bedtools/2.29.2
module load htslib/1.9

# VCF for 8 mutated lines
VCF_SNPs="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/mut8_snps_biallelic.callable.vcf.gz"
#VCF_INDELs="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/mut8_indels_noComplex.callable.vcf.gz"
VCF_INDELs_B="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/mut8_indels_biallelic.callable.vcf.gz"
VCF_INDELs_M="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/mut8_indels_multiallelic.callable.vcf.gz"

OUT_PREFIX="mut8"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered"

# BED file containing sites that differ between 10x morex-sample2 and Morex reference
# SNPs and Indels called in the phased variants VCF
PHV_MOREX_DIFFS_SNPs="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/diffs_from_ref/morex-sample2_phased_variants-snps.callable.biallelic.diffs_from_ref.bed"
PHV_MOREX_DIFFS_INDELS_1bp="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/diffs_from_ref/morex-sample2_phased_variants-indels.callable.biallelic.1bp.diffs_from_ref.bed"
# Use VCF instead of bed where size of variant is correct because in the VCF it's represented as a single position,
# so the range in a bed removes too many variants given how the indel is represented
PHV_MOREX_DIFFS_INDELS_gt1bp="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/diffs_from_ref/morex-sample2_phased_variants-indels.callable.biallelic.gt1bp.vcf.gz"

# Known SNPs
SNPs_BOPA="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k/bopa_idt95_noRescuedSNPs_partsRef.vcf"
SNPs_9k="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k/9k_idt95_noRescuedSNPs_partsRef.vcf"
SNPs_50k="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k/50k_idt95_noRescuedSNPs_partsRef.vcf"

#-------------------
mkdir -p ${OUT_DIR}

function check_filter_stringency() {
    local vcf_file="$1"
    local num_sites="$2"
    if [[ ${num_sites} == 0 ]]; then
        echo "No sites left after filtering. Try using less stringent criteria. File with no sites is: ${vcf_file}. Exiting..." >&2
        exit 8 # If not sites left, error out with message
    fi
}

export -f check_filter_stringency

function count_sites() {
    local vcf="$1"
    local log_file="$2"
    # Get the number of sites left after filtering
    if [[ "${vcf}" == *"gz"* ]]; then
        # We are working with gzipped VCF file
        num_sites=$(zgrep -v "#" ${vcf} | wc -l)
    else
        # We are working with uncompressed vcf
        num_sites=$(grep -v "#" ${vcf} | wc -l)
    fi
    # Append the number of sites remaining to file
    printf "${vcf}\t${num_sites}\n" >> ${log_file}
    check_filter_stringency ${vcf} ${num_sites}
}

export -f count_sites

# Exclude sites that differ between 10x Morex and Morex reference for mutated line samples
# snps
bedtools intersect -v -header -a ${VCF_SNPs} -b ${PHV_MOREX_DIFFS_SNPs} ${SNPs_BOPA} ${SNPs_9k} ${SNPs_50k} | bgzip > ${OUT_DIR}/${OUT_PREFIX}_biallelic.callable.SNPs.noMorexDiffs.vcf.gz
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}_biallelic.callable.SNPs.noMorexDiffs.vcf.gz

# indels
# bedtools intersect -v -header -a ${VCF_INDELs} -b ${PHV_MOREX_DIFFS_INDELS} | bgzip > ${OUT_DIR}/${OUT_PREFIX}_biallelic.callable.INDELs.noMorexDiffs.vcf.gz
# tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}_biallelic.callable.INDELs.noMorexDiffs.vcf.gz
# Indels - biallelic
bedtools intersect -v -header -a ${VCF_INDELs_B} -b ${PHV_MOREX_DIFFS_INDELS_1bp} ${PHV_MOREX_DIFFS_INDELS_gt1bp} | bgzip > ${OUT_DIR}/${OUT_PREFIX}_biallelic.callable.INDELs.noMorexDiffs.vcf.gz
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}_biallelic.callable.INDELs.noMorexDiffs.vcf.gz
# Indels - multiallelic
bedtools intersect -v -header -a ${VCF_INDELs_M} -b ${PHV_MOREX_DIFFS_INDELS_1bp} ${PHV_MOREX_DIFFS_INDELS_gt1bp} | bgzip > ${OUT_DIR}/${OUT_PREFIX}_multiallelic.callable.INDELs.noMorexDiffs.vcf.gz
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}_multiallelic.callable.INDELs.noMorexDiffs.vcf.gz

# Get the number of sites
count_sites ${VCF_SNPs} ${OUT_DIR}/${OUT_PREFIX}_num_sites.log
count_sites ${VCF_INDELs} ${OUT_DIR}/${OUT_PREFIX}_num_sites.log
count_sites ${OUT_DIR}/${OUT_PREFIX}_biallelic.callable.SNPs.noMorexDiffs.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_num_sites.log
#count_sites ${OUT_DIR}/${OUT_PREFIX}_biallelic.callable.INDELs.noMorexDiffs.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_num_sites.log
count_sites ${OUT_DIR}/${OUT_PREFIX}_biallelic.callable.INDELs.noMorexDiffs.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_num_sites.log
count_sites ${OUT_DIR}/${OUT_PREFIX}_multiallelic.callable.INDELs.noMorexDiffs.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_num_sites.log
