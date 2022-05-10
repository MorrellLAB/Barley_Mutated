#!/bin/bash

set -e
set -o pipefail

# Add sample name for VCF so that we can merge multiple samples later on

# Dependencies
module load bcftools/1.10.2
module load htslib/1.9

# User provided input arguments
VCF="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3/Morex_ont_partsRefv3_90.vcf.gz"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3"
NEW_NAME="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3/new_name_ont.txt"

#-------------------
# Make output directories
mkdir -p ${OUT_DIR}

# Rename sample name in VCF and
OUT_PREFIX=$(basename ${VCF} .vcf.gz)
bgzip -dc ${VCF} | bcftools reheader -s ${NEW_NAME} > ${OUT_DIR}/${OUT_PREFIX}_renamed.vcf
