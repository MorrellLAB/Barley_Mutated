#!/bin/bash
#PBS -l mem=12gb,nodes=1:ppn=8,walltime=01:00:00
#PBS -m abe
#PBS -M liux1299@umn.edu
#PBS -q mesabi

set -e
set -o pipefail

# Dependencies
module load bcftools/1.9

# User provided input arguments
# List of VCF files to merge
VCF1=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/M01-3-3/M01-3-3_concat.vcf.gz
VCF2=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/M20-2-2/M20-2-2_concat.vcf.gz
VCF3=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/M29-2-2/M29-2-2_concat.vcf.gz
OUT_DIR=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated
PREFIX="mut_3_lines"

# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR}

# Merge VCF files
bcftools merge ${VCF1} ${VCF2} ${VCF3} -o ${OUT_DIR}/${PREFIX}.vcf.gz -O z --threads 12

# Index concatenated files
bcftools index ${OUT_DIR}/${PREFIX}.vcf.gz -t -o ${OUT_DIR}/${PREFIX}.vcf.gz.tbi
