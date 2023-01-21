#!/bin/bash

set -e
set -o pipefail

# Dependencies
module load bcftools/1.10.2
module load htslib/1.9

# User provided input arguments
M01_VCF="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/M01_ont_partsRefv3/M01_ont_partsRefv3_90_renamed.nochrUn.vcf"
M20_VCF="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/M20_ont_partsRefv3/M20_ont_partsRefv3_90_renamed.nochrUn.vcf"
M29_VCF="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/M29_ont_partsRefv3/M29_ont_partsRefv3_90_renamed.nochrUn.vcf"

OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley"
OUT_PREFIX="mut_ont"

#-------------------
# Compress VCFs before merge
bgzip ${M01_VCF}
bgzip ${M20_VCF}
bgzip ${M29_VCF}
# Index
tabix -p vcf ${M01_VCF}.gz
tabix -p vcf ${M20_VCF}.gz
tabix -p vcf ${M29_VCF}.gz

bcftools merge ${M01_VCF}.gz ${M20_VCF}.gz ${M29_VCF}.gz -m both -O z -o ${OUT_DIR}/${OUT_PREFIX}.vcf.gz
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}.vcf.gz
