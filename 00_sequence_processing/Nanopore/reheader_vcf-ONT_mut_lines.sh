#!/bin/bash

set -e
set -o pipefail

# Add sample name for VCF so that we can merge multiple samples later on

# Dependencies
module load bcftools/1.10.2
module load htslib/1.9

# User provided input arguments
VCF1="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/M01_ont_partsRefv3/M01_ont_partsRefv3_90.vcf.gz"
VCF2="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/M20_ont_partsRefv3/M20_ont_partsRefv3_90.vcf.gz"
VCF3="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/M29_ont_partsRefv3/M29_ont_partsRefv3_90.vcf.gz"
# New names (in same order as above)
NEW_NAME1="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/M01_ont_partsRefv3/new_sample_name.txt"
NEW_NAME2="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/M20_ont_partsRefv3/new_sample_name.txt"
NEW_NAME3="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/M29_ont_partsRefv3/new_sample_name.txt"

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
    bgzip -dc ${vcf} | bcftools reheader -s ${new_name} | bcftools view --targets "^chrUn" > ${out_dir}/${out_prefix}_renamed.nochrUn.vcf
}

export -f rename_vcf_sample

# Rename samples
rename_vcf_sample ${VCF1} ${NEW_NAME1}
rename_vcf_sample ${VCF2} ${NEW_NAME2}
rename_vcf_sample ${VCF3} ${NEW_NAME3}
