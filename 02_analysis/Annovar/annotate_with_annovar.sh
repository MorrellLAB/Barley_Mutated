#!/bin/bash

set -e
set -o pipefail

# This script converst VCF to Annovar's input file format

# Dependencies
module load perl/5.26.1
# Export paths to directories containing annovar scripts so scripts are calleble from anywhere without specifying the path
export PATH=${PATH}:/panfs/jay/groups/9/morrellp/shared/Software/annovar
export PATH=${PATH}:/panfs/jay/groups/9/morrellp/shared/Software/annovar_conversion_tools

# User provided input arguments
annovar_input="$1"
# Full filepath to output directory
#   This directory should also contain the database fasta file
out_dir="$2"
# This must match the *refGene.txt file
#   Example: HV_Morex_v2_HC_refGene.txt should have build version HV_Morex_v2_HC
build_version="$3"

#-----------------
# Check if out directory exists, if not make it
mkdir -p ${out_dir}

cd ${out_dir}

annotate_variation.pl --geneanno --dbtype refGene --buildver ${build_version} ${annovar_input} ${out_dir}/
