#!/bin/bash

set -e
set -o pipefail

# Pull out rare variants (variants with allele count of two or less) from hybrid lines VCF

# Dependencies
module load bcftools/1.10.2
module load htslib/1.9
module load parallel/20210822

# User provided input arguments
VCF_SNPs="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/hybrid13_snps_biallelic.callable.vcf.gz"
VCF_INDELs="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/hybrid13_indels_biallelic.callable.vcf.gz"

OUT_PREFIX="hybrid13"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs"

SAMPLE_LIST="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/hybrid13_sample_names.txt"

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

function split_by_sample() {
    local vcf="$1"
    local curr_sample="$2"
    local out_dir="$3"
    # Prepare subdirectory
    mkdir -p ${out_dir}
    # Prep output file prefix
    if [[ "${vcf}" == *".gz"* ]]; then
        prefix=$(basename ${vcf} .vcf.gz)
    else
        prefix=$(basename ${vcf} .vcf)
    fi
    # Pull out current sample
    # Remove sites where there are missing genotypes or ref-ref genotypes
    bcftools view --samples "${curr_sample}" ${vcf} | bcftools view -e "GT='mis' | GT='RR'" -O z -o ${out_dir}/${curr_sample}_${prefix}.vcf.gz
    tabix -p vcf ${out_dir}/${curr_sample}_${prefix}.vcf.gz
}

export -f split_by_sample

# Get the number of sites for starting VCF and append to file
count_sites ${VCF_SNPs} ${OUT_DIR}/${OUT_PREFIX}_SNPs_num_sites.log
count_sites ${VCF_INDELs} ${OUT_DIR}/${OUT_PREFIX}_INDELs_num_sites.log

# Separate rare vs common variants
# SNPs
bcftools view -i "INFO/AC<=2" ${VCF_SNPs} -O z -o ${OUT_DIR}/${OUT_PREFIX}.SNPs.rare.vcf.gz
bcftools view -i "INFO/AC>=3" ${VCF_SNPs} -O z -o ${OUT_DIR}/${OUT_PREFIX}.SNPs.common.vcf.gz
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}.SNPs.rare.vcf.gz
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}.SNPs.common.vcf.gz
# Count number of sites
count_sites ${OUT_DIR}/${OUT_PREFIX}.SNPs.rare.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_SNPs_num_sites.log
count_sites ${OUT_DIR}/${OUT_PREFIX}.SNPs.common.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_SNPs_num_sites.log

# Indels
bcftools view -i "INFO/AC<=2" ${VCF_INDELs} -O z -o ${OUT_DIR}/${OUT_PREFIX}.INDELs.rare.vcf.gz
bcftools view -i "INFO/AC>=3" ${VCF_INDELs} -O z -o ${OUT_DIR}/${OUT_PREFIX}.INDELs.common.vcf.gz
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}.INDELs.rare.vcf.gz
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}.INDELs.common.vcf.gz
# Count number of sites
count_sites ${OUT_DIR}/${OUT_PREFIX}.INDELs.rare.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_INDELs_num_sites.log
count_sites ${OUT_DIR}/${OUT_PREFIX}.INDELs.common.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_INDELs_num_sites.log

# Pull out homozygous (alt-alt hom) sites for rare variants
bcftools view -i 'GT="AA"' ${OUT_DIR}/${OUT_PREFIX}.SNPs.rare.vcf.gz -O z -o ${OUT_DIR}/${OUT_PREFIX}.SNPs.rare.HOM.vcf.gz
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}.SNPs.rare.HOM.vcf.gz

bcftools view -i 'GT="AA"' ${OUT_DIR}/${OUT_PREFIX}.INDELs.rare.vcf.gz -O z -o ${OUT_DIR}/${OUT_PREFIX}.INDELs.rare.HOM.vcf.gz
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}.INDELs.rare.HOM.vcf.gz

# Pull out heterozygous sites only
bcftools view -i 'GT="het"' ${OUT_DIR}/${OUT_PREFIX}.SNPs.rare.vcf.gz -O z -o ${OUT_DIR}/${OUT_PREFIX}.SNPs.rare.HET.vcf.gz
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}.SNPs.rare.HET.vcf.gz

bcftools view -i 'GT="het"' ${OUT_DIR}/${OUT_PREFIX}.INDELs.rare.vcf.gz -O z -o ${OUT_DIR}/${OUT_PREFIX}.INDELs.rare.HET.vcf.gz
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}.INDELs.rare.HET.vcf.gz

# Check number of sites
count_sites ${OUT_DIR}/${OUT_PREFIX}.SNPs.rare.HOM.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_SNPs_num_sites.log
count_sites ${OUT_DIR}/${OUT_PREFIX}.SNPs.rare.HET.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_SNPs_num_sites.log
count_sites ${OUT_DIR}/${OUT_PREFIX}.INDELs.rare.HOM.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_INDELs_num_sites.log
count_sites ${OUT_DIR}/${OUT_PREFIX}.INDELs.rare.HET.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_INDELs_num_sites.log

# Split by sample
# Some downstream analyses work better when each VCF only has one sample without
#   a bunch of extra missing genotypes
parallel --verbose split_by_sample ${VCF_SNPs} {} "${OUT_DIR}/split_by_sample_SNPs" :::: ${SAMPLE_LIST}
parallel --verbose split_by_sample ${OUT_DIR}/${OUT_PREFIX}.SNPs.rare.vcf.gz {} "${OUT_DIR}/split_by_sample_SNPs_rare" :::: ${SAMPLE_LIST}
parallel --verbose split_by_sample ${OUT_DIR}/${OUT_PREFIX}.SNPs.common.vcf.gz {} "${OUT_DIR}/split_by_sample_SNPs_common" :::: ${SAMPLE_LIST}

parallel --verbose split_by_sample ${VCF_INDELs} {} "${OUT_DIR}/split_by_sample_INDELs" :::: ${SAMPLE_LIST}
parallel --verbose split_by_sample ${OUT_DIR}/${OUT_PREFIX}.INDELs.rare.vcf.gz {} "${OUT_DIR}/split_by_sample_INDELs_rare" :::: ${SAMPLE_LIST}
parallel --verbose split_by_sample ${OUT_DIR}/${OUT_PREFIX}.INDELs.common.vcf.gz {} "${OUT_DIR}/split_by_sample_INDELs_common" :::: ${SAMPLE_LIST}
