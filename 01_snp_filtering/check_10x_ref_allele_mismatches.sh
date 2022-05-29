#!/bin/bash

set -e
set -o pipefail

# User provided input arguments
REF="$1"
DEL_VCF="$2"
LARGE_SV_VCF="$3"
PHASED_VAR_VCF="$4"

#-------------------
# Get current sample name
sample_name=$(basename ${DEL_VCF} _dels.vcf.gz)

# DEL
bcftools norm -c w -f ${REF} -o temp_ref_check_${sample_name}_dels.vcf -O v ${DEL_VCF} >& temp_ref_check_warn_${sample_name}_dels.log
# Check unique warnings
echo "Unique warnings for DELs:"
cut -f 1 temp_ref_check_warn_${sample_name}_dels.log | sort -u

# Large SVs
bcftools norm -c w -f ${REF} -o temp_ref_check_${sample_name}_large_svs.vcf -O v ${LARGE_SV_VCF} >& temp_ref_check_warn_${sample_name}_large_svs.log
# Check unique warnings
printf "\n"
echo "Unique warnings for Large SVs:"
cut -f 1 temp_ref_check_warn_${sample_name}_large_svs.log | sort -u

# Phased variants
bcftools norm -c w -f ${REF} -o temp_ref_check_${sample_name}_phased_variants.vcf -O v ${PHASED_VAR_VCF} >& temp_ref_check_warn_${sample_name}_phased_variants.log
# Check unique warnings
printf "\n"
echo "Unique warnings for Phased Variants:"
cut -f 1 temp_ref_check_warn_${sample_name}_phased_variants.log | sort -u
