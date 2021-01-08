#!/bin/bash

set -e
set -o pipefail

# This scripts runs BAD_Mutations VeP_to_Subs.py script to convert the
#   VeP reports into a format for BAD_Mutations to compile predictions

module load parallel
module load python3/3.6.3_anaconda5.0.1

# Activate conda environment
source activate /home/morrellp/liux1299/.conda/envs/bad_mutations

# User defined input arguments
VEP_REPORT_LIST=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/HC_parts_gff/for_bad_mut/vep_report_list.txt
# Full path to out directory
OUT_DIR=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/vep_to_subs
# Path to VeP_to_Subs.py script directory
SCRIPT_DIR=~/Software/BAD_Mutations/Supporting

#------------------------------
# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR}

# Store list of VeP reports in array
VEP_ARR=($(cat ${VEP_REPORT_LIST}))

function vep_to_subs() {
    local vep_report="$1"
    local out_dir="$2"
    local script_dir="$3"
    out_prefix=$(basename ${vep_report} .txt.gz)
    # Make subdirectory for each report
    mkdir -p ${out_dir}/per_transcript_subs_${out_prefix}
    # Run VeP to subs
    python ${script_dir}/VeP_to_Subs.py \
        ${vep_report} \
        ${out_dir}/per_transcript_subs_${out_prefix}/${out_prefix}_long_subs.txt \
        ${out_dir}/per_transcript_subs_${out_prefix}
}

export -f vep_to_subs

parallel vep_to_subs {} ${OUT_DIR} ${SCRIPT_DIR} ::: ${VEP_ARR[@]}
