#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=12gb
#SBATCH --tmp=10gb
#SBATCH -t 03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# This script pulls out invariant sites only (which should be filtered separately)

# Dependencies
module load bcftools/1.10.2
module load htslib/1.9

# User provided input arguments
VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/Genotype_GVCFs/morex-sample2.vcf.gz"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/Genotype_GVCFs"
OUT_PREFIX="morex-sample2"

#-----------------
mkdir -p ${OUT_DIR}

# Select invariant sites only
bcftools view -i 'GT[*]="RR"' ${VCF} --threads 12 -O z -o ${OUT_DIR}/${OUT_PREFIX}.invariant_sites.vcf.gz
# Index vcf
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}.invariant_sites.vcf.gz
