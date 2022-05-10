#!/bin/bash

set -e
set -o pipefail

# Add sample name for VCF so that we can merge multiple samples later on

# Dependencies
module load bcftools/1.10.2
module load htslib/1.9

# User provided input arguments
VCF1="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3/Morex_ont_partsRefv3_90.vcf.gz"
VCF2=""
VCF3=""
NEW_NAME1="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3/new_name_ont.txt"
NEW_NAME2=""

#-------------------
# Make output directories
mkdir -p ${OUT_DIR}

# Rename sample name in VCF and
function rename_vcf_sample() {
    local vcf="$1"
    local new_name="$2"
    out_prefix=$(basename ${vcf} .vcf.gz)
    # Output in same directory as current VCF
    out_dir=$(dirname ${vcf})
    bgzip -dc ${vcf} | bcftools reheader -s ${new_name} > ${out_dir}/${out_prefix}_renamed.vcf
}

export -f rename_vcf_sample


