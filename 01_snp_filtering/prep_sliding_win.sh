#!/bin/bash

set -e
set -o pipefail

# Dependencies
module load bedtools/2.29.2

# User provided input arguments
win_size_bp="3000"
step_size_bp="1000"
ref_fai="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta.fai"
out_dir="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/check_heterogeneity"

#-----------------
mkdir -p ${out_dir}

# Prepare sliding windows
bedtools makewindows -g ${ref_fai} -w ${win_size_bp} -s ${step_size_bp} | grep -v "chrUn" > ${out_dir}/sliding_win_${win_size_bp}_overlap_${step_size_bp}.bed

# Reformat regions
sed 's/\t/:/' ${out_dir}/sliding_win_${win_size_bp}_overlap_${step_size_bp}.bed | sed 's/\t/-/' > ${out_dir}/sliding_win_${win_size_bp}_overlap_${step_size_bp}_reformatted.txt
