#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=12gb
#SBATCH --tmp=8gb
#SBATCH -t 02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# This script concatenates a list of split by chromosome VCF files

# Dependencies
module load bcftools/1.10.2
module load htslib/1.9

# User provided input arguments
VCF_LIST="/scratch.global/liux1299/temp_morex-sample2/vcf_split_chr_parts_list.txt"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/pixy_input"
OUT_PREFIX="morex-sample2"

#-----------------
mkdir -p ${OUT_DIR}

# Concatenate split vcf files
bcftools concat --file-list ${VCF_LIST} --threads 12 -O z -o ${OUT_DIR}/${OUT_PREFIX}.vcf.gz
# Index vcf
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}.vcf.gz
