#!/bin/bash

set -e
set -o pipefail

# Add sample name for VCF so that we can merge multiple samples later on

# Dependencies
module load bcftools/1.10.2
module load htslib/1.9

# User provided input arguments
VCF="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_85x_ont/Morex_85X_sorted.vcf"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_85x_ont"
NEW_NAME="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_85x_ont/new_name_ont.txt"

#-------------------
# Make output directories
mkdir -p ${OUT_DIR}

# Rename sample name in VCF and
OUT_PREFIX=$(basename ${VCF} .vcf)
bcftools reheader -s ${NEW_NAME} ${VCF} > ${OUT_DIR}/${OUT_PREFIX}_renamed.vcf
