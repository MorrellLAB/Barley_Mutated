#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=10gb
#SBATCH --tmp=2gb
#SBATCH -t 48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o bad_mut_predict.sh.%A_%a.out
#SBATCH -e bad_mut_predict.sh.%A_%a.err

set -e
set -o pipefail

# This script stores all filepaths and calls on the script bad_mut_predict.sh

# Dependencies
module load parallel
module load python3/3.8.3_anaconda2020.07_mamba

# Activate conda environment
source activate /home/morrellp/liux1299/.conda/envs/bad_mutations

# User provided input arguments
# Full path to a list of lists (to utilize GNU parallel and job arrays)
# List of lists here should only include FASTA files that had an alignment
#	(i.e., .fa and .tree files were generated) in the previous align step
FASTA_LIST_OF_LISTS="/panfs/jay/groups/9/morrellp/shared/Projects/WBDC_inversions/bad_mutations/final_lists/all_cds_hvulgare_list_of_lists.txt"

# Full path to the config file
#CONFIG_FILE="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/config.txt"
CONFIG_FILE="/panfs/jay/groups/9/morrellp/gfrascar/bad_mutations_scripts/config.txt"

# Full path to a list of MSA_Output directories that contain *.fa and *.tree files
#MSA_DIR_LIST=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/MSA_output/all_msa_output_dir_list.txt
MSA_DIR_LIST="/panfs/jay/groups/9/morrellp/shared/Projects/WBDC_inversions/bad_mutations/results/MSA_output/MSA_out_dir_list_of_dir.txt"

# Full path to per transcript substitutions .subs files
#	This output is from the VeP_to_Subs.py supporting script
SUBS_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/vep_to_subs-mut/per_transcript_subs-mut_SNPs_private"

# Full path to a list of primary transcripts, one per line
#PRIMARY_TRANSCRIPTS=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/hvulgare_primary_transcripts_only.txt
PRIMARY_TRANSCRIPTS="/panfs/jay/groups/9/morrellp/shared/Projects/WBDC_inversions/bad_mutations/results/phytozome13_download_V3_primary_transcript/hvulgare_primary_transcripts_only.txt"

# Sample name will be used as a prefix for outputs
SAMPLE_NAME="mut_lines"

# Full path to output directory
OUT_DIR="/scratch.global/liux1299/bad_mutations/predict_output_${SAMPLE_NAME}"

# Full path to the BAD_Mutations.py script
BAD_MUT_SCRIPT="/panfs/jay/groups/9/morrellp/liux1299/Software/BAD_Mutations/BAD_Mutations.py"

# Full path to where we want to store the log files output from parallel
LOG_FILE_DIR="${OUT_DIR}/all_parallel_log_files"

# Script that stores predict_sub function
PREDICT_SCRIPT="/panfs/jay/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/02_analysis/bad_mutations/bad_mut_predict.sh"

#------------------------------
# Run predict script
${PREDICT_SCRIPT} ${FASTA_LIST_OF_LISTS} \
	${CONFIG_FILE} \
	${MSA_DIR_LIST} \
	${SUBS_DIR} \
	${PRIMARY_TRANSCRIPTS} \
	${SAMPLE_NAME} \
	${OUT_DIR} \
	${BAD_MUT_SCRIPT} \
	${LOG_FILE_DIR}
