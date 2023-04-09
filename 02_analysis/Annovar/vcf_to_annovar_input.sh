#!/bin/bash

set -e
set -o pipefail

# This script convert VCF to Annovar's input file format

# Dependencies
module load perl/5.26.1
# Export paths to directories containing annovar scripts so scripts are calleble from anywhere without specifying the path
export PATH=${PATH}:/panfs/jay/groups/9/morrellp/shared/Software/annovar
export PATH=${PATH}:/panfs/jay/groups/9/morrellp/shared/Software/annovar_conversion_tools

# User provided input arguments
vcf_file="$1"
# Full filepath to output directory
out_dir="$2"

#-----------------
# Generate output file pfrefix
out_prefix=$(basename ${vcf_file} .vcf.gz)

convert2annovar.pl --format vcf4 --allsample --withfreq --includeinfo --outfile ${out_dir}/${out_prefix}_annovar_input.txt ${vcf_file}
