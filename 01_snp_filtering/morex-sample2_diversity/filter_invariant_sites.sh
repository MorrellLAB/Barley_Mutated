#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=12gb
#SBATCH --tmp=10gb
#SBATCH -t 08:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# This script filters vcf given cutoffs

# Dependencies
module load bcftools/1.10.2
module load htslib/1.9

# User provided input arguments
VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/Genotype_GVCFs/morex-sample2.invariant_sites.vcf.gz"
OUT_PREFIX="morex-sample2"
# Output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/Genotype_GVCFs/filtered"
# Temporary directory
TEMP_DIR="/scratch.global/liux1299/temp_morex-sample2"

# Mininimum per sample DP threshold
MIN_DP="5"
# Quality cutoff
QUAL_CUTOFF="30"
# Downsampled chr1H, 6th percentile is 4
RGQ_CUTOFF="5"

#-----------------
mkdir -p ${OUT_DIR} ${TEMP_DIR}

# Select invariant sites only
# Remove sites where REF column is missing
bcftools filter -e "ALT='<NON_REF>' | FMT/DP<${MIN_DP} | FMT/RGQ<${RGQ_CUTOFF} | QUAL<${QUAL_CUTOFF}" ${VCF} -O z -o ${OUT_DIR}/${OUT_PREFIX}.invariant_sites.filtDP-RGQ-QUAL.vcf.gz
# Index vcf
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}.invariant_sites.filtDP-RGQ-QUAL.vcf.gz

# Pull out sites where REF column is missing into temporary file
# For exploration purposes and to check some sites
bcftools view -i "ALT='<NON_REF>'" ${VCF} -O v -o ${TEMP_DIR}/${OUT_PREFIX}.invariant_sites.nonREF.vcf
