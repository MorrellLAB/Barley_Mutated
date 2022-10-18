#!/bin/bash

# This script checks the number of segregating sites in a provided region

# Usage: ./check_seg_sites.sh <QUERY_REGION> <SANGER_VCF> <OUTFILE_PREFIX> {VCF} {OUT_DIR}
# Where:
#   <QUERY_REGION> is our query Sanger region formatted for bcftools --regions flag, chr:start-stop
#   <SANGER_VCF> is the path to the Sanger region VCF file
#   <OUTFILE_PREFIX> should be the gene/locus name, example: Adh1
#   {VCF} has a default value that is used if no arguments were passed. This is the path to the VCF with SNP calls
#   {OUT_DIR} has a default value that is used if no arguments were passed. This is the path to our output directory
# Note: <required_argument> {default_values}

# Dependencies
module load bcftools/1.10.2
module load bedtools/2.29.2

# User provided input arguments
# Prepare default values to use if not specified
DEFAULT_VCF="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/mut8_and_hybrid_barley_snps_biallelic.noRepeatOverlap.noRefNs.vcf.gz"
DEFAULT_OUTDIR="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/check_seg_sites_and_common_snps"

# Our query Sanger region formatted for bcftools --regions flag, chr:start-stop
REGION="$1"
SANGER_VCF="$2"
# Output file prefix should be the gene/locus name
OUT_PREFIX="$3"
# Use default value if no arguments were passed, otherwise use provided input argument
# VCF with SNP calls
VCF="${4:-$DEFAULT_VCF}"
# Output directory
OUT_DIR="${5:-$DEFAULT_OUTDIR}"

#------------------------
# Make output directory
mkdir -p ${OUT_DIR}

# Prepare VCF basename
if [[ ${VCF} == *".vcf.gz" ]]; then
    VCF_BN=$(basename ${VCF} .vcf.gz)
elif [[ ${VCF} == *".vcf" ]]; then
    VCF_BN=$(basename ${VCF} .vcf)
fi

# Save segregating sites to file
bcftools view --regions "${REGION}" ${VCF} -O v -o ${OUT_DIR}/${OUT_PREFIX}-${VCF_BN}.vcf

# Count the number of SNPs in query region
num_sites_in_region=$(grep -v "#" ${OUT_DIR}/${OUT_PREFIX}-${VCF_BN}.vcf | wc -l)

if [[ ${num_sites_in_region} == "0" ]]; then
    # No sites from VCF in query Sanger region
    echo "No sites found in region ${REGION} in VCF: ${VCF}"
    # Delete empty VCF
    rm ${OUT_DIR}/${OUT_PREFIX}-${VCF_BN}.vcf
else
    echo "${num_sites_in_region} sites found in region ${REGION} in VCF: ${VCF}"
    # Pull intersecting Sanger sites
    bedtools intersect -header -wa -a ${SANGER_VCF} -b ${OUT_DIR}/${OUT_PREFIX}-${VCF_BN}.vcf > ${OUT_DIR}/${OUT_PREFIX}-sanger_x_${VCF_BN}.vcf
    # Check counts in sanger intersecting query region vcf
    num_sanger_x_vcf=$(grep -v "#" ${OUT_DIR}/${OUT_PREFIX}-sanger_x_${VCF_BN}.vcf | wc -l)
    echo "${num_sanger_x_vcf} sites found in ${SANGER_VCF}"
fi
