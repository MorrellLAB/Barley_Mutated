#!/bin/bash

set -e
set -o pipefail

# Dependencies
module load bcftools/1.10.2
module load htslib/1.9
module load parallel/20210822

# Split into one sample per file, easier for downstream reformatting for plotting
VCF_10x_dels="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/de_novo_larger_svs/mut_3_lines_dels_merged.callable.noRefDiffs.private.supports.vcf"
# Output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/de_novo_larger_svs"
# Sample names must match VCF sample names exactly
SAMPLE_ARR=("M01-3-3" "M20-2-2" "M29-2-2")

#------------------
mkdir -p ${OUT_DIR}

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

parallel --verbose split_by_sample ${VCF_10x_dels} {} ${OUT_DIR} ::: ${SAMPLE_ARR[@]}

# Rename files so they are shorter
rename -v "_mut_3_lines_dels_merged" "" *_mut_3_lines_dels_merged.callable.noRefDiffs.private.supports.vcf.gz*
