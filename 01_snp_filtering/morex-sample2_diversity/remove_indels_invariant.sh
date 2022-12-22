#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=12gb
#SBATCH --tmp=8gb
#SBATCH -t 06:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# This script concatenates a list of split by chromosome VCF files

# Dependencies
module load perl/5.26.1
module load vcflib_ML/1.0.0_rc2

# User provided input arguments
VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/Genotype_GVCFs/filtered/morex-sample2.invariant_sites.filtDP-RGQ-QUAL.vcf.gz"
# Temporary directory
TEMP_DIR="/scratch.global/liux1299/temp_morex-sample2"

#-----------------
mkdir -p ${TEMP_DIR}

# Get output file basename and output directory
BN=$(basename ${VCF} .vcf.gz)

# Exclude indels
zcat ${VCF} | vcfnoindels > ${TEMP_DIR}/${BN}.snps.vcf

# Pull out indels only (for checking)
zcat ${VCF} | vcfindels > ${TEMP_DIR}/${BN}.indels.vcf
