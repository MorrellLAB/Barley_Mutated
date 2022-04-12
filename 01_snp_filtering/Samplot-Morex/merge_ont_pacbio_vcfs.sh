#!/bin/bash

# This script serves as a log of commands run

set -e
set -o pipefail

# Dependencies
module load bcftools/1.9

# User provided input arguments
# List of VCF files to merge
VCF_ONT="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3/filtered/morex_ont_noHomRef.vcf.gz"
VCF_PACBIO="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/pacbio_morex/pacbio_morex_v3/filtered/morex_pacbio_noHomRef.vcf.gz"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration_morex_v3"
PREFIX="morex_ont_pacbio_noHomRef"

#---------------------
# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR}

# Merge VCF files
bcftools merge ${VCF_ONT} ${VCF_PACBIO} -o ${OUT_DIR}/${PREFIX}.vcf.gz -O z --threads 6
# Index VCF
tabix -p vcf ${OUT_DIR}/${PREFIX}.vcf.gz
