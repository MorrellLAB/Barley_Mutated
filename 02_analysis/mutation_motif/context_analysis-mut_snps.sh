#!/bin/bash

set -e
set -o pipefail

# Dependencies
module load parallel/20210822
module load python3/3.8.3_anaconda2020.07_mamba
# Conda environment for mutation motif software
source activate /panfs/jay/groups/9/morrellp/liux1299/.conda/envs/mut_motif_env
# Load shared functions
source ~/GitHub/Barley_Mutated/02_analysis/mutation_motif/context_analysis_functions.sh

# User provided input arguments
# Combined counts files
COMBINED_COUNTS_FILE="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/mutation_motif/mut8_and_3mut10xGenomics.SNPs.private/counts_tables/combined_counts_mut8_and_3mut10xGenomics.SNPs.private.winFlank100bp.txt"
# List of counts separated by mutation direction
COUNTS_LIST="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/mutation_motif/mut8_and_3mut10xGenomics.SNPs.private/counts_tables/separate_counts_file_list.txt"
# Parent output directory (subdirectories will automatically be created)
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/mutation_motif/mut8_and_3mut10xGenomics.SNPs.private"
# Used to rename output files so they are easier to collate later on
# Example: the output file "1.pdf" in the AtoT subdirectory would get renamed to "1_AtoT_${PLOT_FILE_BN}.pdf"
PLOT_FILE_BN="mut.SNPs.private"
# Must match portion of counts filename exactly. Will be used to autogenerate output subdirectories
#    Example filename: AtoC_mut8_and_3mut10xGenomics.SNPs.private.winFlank100bp.txt
#    COUNTS_BN to use: mut8_and_3mut10xGenomics.SNPs.private.winFlank100bp
COUNTS_BN="mut8_and_3mut10xGenomics.SNPs.private.winFlank100bp"

#------------------
# Make output directory and subdirectories
mkdir -p ${OUT_DIR} \
    ${OUT_DIR}/simple_analysis \
    ${OUT_DIR}/strand_symmetry_analysis \
    ${OUT_DIR}/strand_symmetry_analysis/strand_symmetry_counts \
    ${OUT_DIR}/full_spectra_analysis

# Match what is in the counts table
DIRECTION_ARR=("AtoC" "AtoG" "AtoT" "CtoA" "CtoG" "CtoT" "GtoA" "GtoC" "GtoT" "TtoA" "TtoC" "TtoG")

# Evaluating the effect of neighbours on mutation
# Simple analysis
# Also generates plots
parallel --verbose simple_analysis {} ${OUT_DIR}/simple_analysis ${COUNTS_BN} ${PLOT_FILE_BN} :::: ${COUNTS_LIST}

# Test for strand symmetry (or asymmetry)
# Also generates plots
parallel --verbose test_strand_symmetry {} ${COMBINED_COUNTS_FILE} "${OUT_DIR}/strand_symmetry_analysis" ${COUNTS_BN} ${PLOT_FILE_BN} ::: ${DIRECTION_ARR[@]}

# Test full spectra
# Also generates plots
full_spectra_analysis ${COMBINED_COUNTS_FILE} ${OUT_DIR}/full_spectra_analysis ${PLOT_FILE_BN}
