#!/bin/bash

set -e
set -o pipefail

# This script renames the file extension .fa to .fasta

# User provided input arguments
ALL_CDS_FASTA_LIST=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/align_lists/all_cds_hvulgare_list.txt
MSA_DIR_LIST=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/MSA_output/all_msa_output_dir_list.txt

#----------------------------
# Rename from .fa to .fasta
for i in $(cat ${ALL_CDS_FASTA_LIST})
do
    # Check fasta file extension is .fasta and NOT .fa
    # If not .fasta, rename
    if [[ ${i} != *".fasta"* ]]; then
        rename -v ".fa" ".fasta" ${i}
    fi
done

# Repeat for MSA output directory fasta files
# Build array of MSA_output subdirectories containing .fa and .tree files
MSA_DIR_ARR=($(cat ${MSA_DIR_LIST}))
for i in ${MSA_DIR_ARR[@]}
do
    rename -v ".fa" ".fasta" ${i}/*.fa
done
