#!/bin/bash

set -e
set -o pipefail

# Run all steps to determine if BAD_Mutations predictions are "Deleterious" vs "Tolerated"

# Dependencies
module load R/4.1.0
module load python3/3.8.3_anaconda2020.07_mamba
# Path to directory containing scripts called below
export PATH=${PATH}:~/GitHub/Barley_Mutated/02_analysis/bad_mutations

# User provided input arguments
# Compiled predictions after final step in BAD_mutations pipeline
BAD_MUT_REPORT="$1"
VEP_REPORT="$2"
OUT_PREFIX="$3"
# Criteria for calling SNP "Deleterious"
MIN_SEQ="$4"
MAX_CONSTRAINT="$5"
P_CUTOFF="$6"
# Note: this should be the number of nonsynonymous SNPs given as input to the BAD_Mutations predict step
#   and NOT the number of lines in the compiled report
NUM_CODONS_TESTED="$7"

#-------------------
# Prepare output directories and file prefix
bad_mut_report_bn=$(basename ${BAD_MUT_REPORT} .txt)
bad_mut_report_basedir=$(dirname ${BAD_MUT_REPORT})

# Merge BAD_Mutations compiled predictions and VeP report
merge_bad_mut_and_vep.R ${BAD_MUT_REPORT} ${VEP_REPORT} ${bad_mut_report_basedir}/${bad_mut_report_bn}_with_VeP.txt

# Annotate SNP as "Deleterious" vs "Tolerated" using heuristic approach
dSNP_BADMutation_heuristics_masked-vep.py ${bad_mut_report_basedir}/${bad_mut_report_bn}_with_VeP.txt ${MIN_SEQ} ${MAX_CONSTRAINT} ${P_CUTOFF} ${NUM_CODONS_TESTED} > ${bad_mut_report_basedir}/${OUT_PREFIX}_deleterious_vs_tolerated.txt
