#!/bin/bash

set -e
set -o pipefail

# Combine SNPs and 1 bp indel mutated samples into one VCF and identify variants
#   private to each sample

# Dependencies
module load bcftools/1.10.2
module load htslib/1.9

# User provided input arguments
# 8 WGS mutated lines VCF SNPs and 1 bp indels
MUT8_VCF_SNPs="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/mut8_biallelic.noRepeatOverlap.noRefNs.SNPs.noMorexDiffs.vcf.gz"
MUT8_VCF_INDELs="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/mut8_biallelic.noRepeatOverlap.noRefNs.INDELs.noMorexDiffs.vcf.gz"

# 10x Genomics phased variants VCF SNPs and 1 bp indels
M01_VCF_SNPS="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/M01-3-3_phased_variants.noRepeatOverlap.noRefNs.SNPs.noMorexDiffs.vcf.gz"
M01_VCF_INDELs="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/M01-3-3_phased_variants.noRepeatOverlap.noRefNs.INDELs.noMorexDiffs.vcf.gz"

M20_VCF_SNPs="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/M20-2-2_phased_variants.noRepeatOverlap.noRefNs.SNPs.noMorexDiffs.vcf.gz"
M20_VCF_INDELs="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/M20-2-2_phased_variants.noRepeatOverlap.noRefNs.INDELs.noMorexDiffs.vcf.gz"

M29_VCF_SNPs="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/M29-2-2_phased_variants.noRepeatOverlap.noRefNs.SNPs.noMorexDiffs.vcf.gz"
M29_VCF_INDELs="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/M29-2-2_phased_variants.noRepeatOverlap.noRefNs.INDELs.noMorexDiffs.vcf.gz"

OUT_PREFIX="mut8_and_3mut10xGenomics"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs"

#----------------
mkdir -p ${OUT_DIR} ${OUT_DIR}/Intermediates

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

# Merge mutated line samples into one VCF
# SNPs
bcftools merge \
    ${MUT8_VCF_SNPs} \
    ${M01_VCF_SNPS} \
    ${M20_VCF_SNPs} \
    ${M29_VCF_SNPs} \
    -O z -o ${OUT_DIR}/Intermediates/${OUT_PREFIX}.SNPs.vcf.gz
tabix -p vcf ${OUT_DIR}/Intermediates/${OUT_PREFIX}.SNPs.vcf.gz

# Indels
bcftools merge \
    ${MUT8_VCF_INDELs} \
    ${M01_VCF_INDELs} \
    ${M20_VCF_INDELs} \
    ${M29_VCF_INDELs} \
    -O z -o ${OUT_DIR}/Intermediates/${OUT_PREFIX}.INDELs.vcf.gz
tabix -p vcf ${OUT_DIR}/Intermediates/${OUT_PREFIX}.INDELs.vcf.gz

# Get the number of sites left for starting VCF and append to file
count_sites ${OUT_DIR}/Intermediates/${OUT_PREFIX}.SNPs.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_SNPs_num_sites.log
count_sites ${OUT_DIR}/Intermediates/${OUT_PREFIX}.INDELs.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_INDELs_num_sites.log

# Exclude non-unique variants (i.e., variants that are present in more than 1 of the mutated lines)
# We want variants private to each sample
# First, create a file where only one sample has an alternate allele
# SNPs
bcftools view -i "COUNT(GT='alt')=1" ${OUT_DIR}/Intermediates/${OUT_PREFIX}.SNPs.vcf.gz -O z -o ${OUT_DIR}/${OUT_PREFIX}.SNPs.private.vcf.gz
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}.SNPs.private.vcf.gz
# Check number of sites
count_sites ${OUT_DIR}/${OUT_PREFIX}.SNPs.private.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_SNPs_num_sites.log

# Indels
bcftools view -i "COUNT(GT='alt')=1" ${OUT_DIR}/Intermediates/${OUT_PREFIX}.INDELs.vcf.gz -O z -o ${OUT_DIR}/${OUT_PREFIX}.INDELs.private.vcf.gz
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}.INDELs.private.vcf.gz
# Check number of sites
count_sites ${OUT_DIR}/${OUT_PREFIX}.INDELs.private.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_INDELs_num_sites.log

# Pull out homozygous (alt-alt hom) sites only
bcftools view -i 'GT="AA"' ${OUT_DIR}/${OUT_PREFIX}.SNPs.private.vcf.gz -O z -o ${OUT_DIR}/${OUT_PREFIX}.SNPs.private.HOM.vcf.gz
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}.SNPs.private.HOM.vcf.gz

bcftools view -i 'GT="AA"' ${OUT_DIR}/${OUT_PREFIX}.INDELs.private.vcf.gz -O z -o ${OUT_DIR}/${OUT_PREFIX}.INDELs.private.HOM.vcf.gz
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}.INDELs.private.HOM.vcf.gz

# Pull out heterozygous sites only
bcftools view -i 'GT="het"' ${OUT_DIR}/${OUT_PREFIX}.SNPs.private.vcf.gz -O z -o ${OUT_DIR}/${OUT_PREFIX}.SNPs.private.HET.vcf.gz
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}.SNPs.private.HET.vcf.gz

bcftools view -i 'GT="het"' ${OUT_DIR}/${OUT_PREFIX}.INDELs.private.vcf.gz -O z -o ${OUT_DIR}/${OUT_PREFIX}.INDELs.private.HET.vcf.gz
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}.INDELs.private.HET.vcf.gz

# Check number of sites
count_sites ${OUT_DIR}/${OUT_PREFIX}.SNPs.private.HOM.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_SNPs_num_sites.log
count_sites ${OUT_DIR}/${OUT_PREFIX}.SNPs.private.HET.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_SNPs_num_sites.log
count_sites ${OUT_DIR}/${OUT_PREFIX}.INDELs.private.HOM.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_INDELs_num_sites.log
count_sites ${OUT_DIR}/${OUT_PREFIX}.INDELs.private.HET.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_INDELs_num_sites.log
