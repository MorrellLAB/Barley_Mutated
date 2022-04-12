#!/bin/bash

set -e
set -o pipefail

# Add sample name for VCF so that we can merge multiple samples later on

# Dependencies
module load bcftools/1.10.2
module load htslib/1.9

# User provided input arguments
VCF="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/pacbio_morex/pacbio_morex_v3/Morex_pacbio_90.vcf.gz"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/pacbio_morex/pacbio_morex_v3"
NEW_NAME="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/pacbio_morex/pacbio_morex_v3/new_name_pacbio.txt"

#-------------------
# Make output directories
mkdir -p ${OUT_DIR}

# Rename sample name in VCF and
OUT_PREFIX=$(basename ${VCF} .vcf.gz)
bgzip -dc ${VCF} | bcftools reheader -s ${NEW_NAME} > ${OUT_DIR}/${OUT_PREFIX}_renamed.vcf
