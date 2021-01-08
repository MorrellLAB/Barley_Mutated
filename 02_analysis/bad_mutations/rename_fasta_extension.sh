#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=2gb
#SBATCH -t 03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o rename_fasta_extension.sh.%A_%a.out
#SBATCH -e rename_fasta_extension.sh.%A_%a.err

set -e
set -o pipefail

# This script renames the file extension .fa to .fasta

# User provided input arguments
# List of all fasta files
ALL_CDS_FASTA_LIST=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/align_lists/all_cds_hvulgare_list.txt
# List of MSA_output subdirectories
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
