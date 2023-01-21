#!/bin/bash

set -e
set -o pipefail

# This script prepares VeP output files then runs BAD_Mutations VeP_to_Subs.py script to convert the
#   VeP reports into a format for BAD_Mutations to compile predictions

module load python3/3.8.3_anaconda2020.07_mamba

# User defined input arguments
VEP_REPORT="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/HC_LC_gff_SNPs_hybrid13_common/hybrid13.SNPs.common.txt"
SUBDIR_BN="hybrid_SNPs_common"
# Full path to out directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/vep_to_subs-hybrid13"
# Path to VeP_to_Subs.py script directory
SCRIPT_DIR="/panfs/jay/groups/9/morrellp/liux1299/Software/BAD_Mutations/Supporting"

#------------------------------
# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR} ${OUT_DIR}/per_transcript_subs-${SUBDIR_BN}

# Get dir of vep report
vep_dir=$(dirname ${VEP_REPORT})

# Prepare input file for BAD_Mutations
# Pull nonsynonymous variants as defined by VeP, which uses
#   Sequence Ontology's definition (http://www.sequenceontology.org/)
if [[ "${VEP_REPORT}" == *"gz"* ]]; then
    out_prefix=$(basename ${VEP_REPORT} .txt.gz)
    # Pull out nonsynonymous variants defined by VeP
    zgrep "#" ${VEP_REPORT} > ${vep_dir}/${out_prefix}.nonsyn.txt
    zgrep 'missense_variant\|start_lost\|stop_gained\|stop_lost' ${VEP_REPORT} | uniq >> ${vep_dir}/${out_prefix}.nonsyn.txt
else
    out_prefix=$(basename ${VEP_REPORT} .txt)
    # Pull out nonsynonymous variants defined by VeP
    grep "#" ${VEP_REPORT} > ${vep_dir}/${out_prefix}.nonsyn.txt
    grep 'missense_variant\|start_lost\|stop_gained\|stop_lost' ${VEP_REPORT} | uniq >> ${vep_dir}/${out_prefix}.nonsyn.txt
fi
# Gzip .txt file
gzip ${vep_dir}/${out_prefix}.nonsyn.txt

# Run VeP to subs
python ${SCRIPT_DIR}/VeP_to_Subs.py \
    "${vep_dir}/${out_prefix}.nonsyn.txt.gz" \
    "${OUT_DIR}/per_transcript_subs-${SUBDIR_BN}/${out_prefix}_long_subs.txt" \
    "${OUT_DIR}/per_transcript_subs-${SUBDIR_BN}"
