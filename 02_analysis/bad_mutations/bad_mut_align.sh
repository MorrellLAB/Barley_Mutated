#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=16gb
#SBATCH --tmp=14gb
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t,max
#SBATCH -o bad_mut_align.sh.%A_%a.out
#SBATCH -e bad_mut_align.sh.%A_%a.err

set -e
set -o pipefail

module load parallel
module load python3/3.6.3_anaconda5.0.1

# Activate conda environment
source activate /home/morrellp/liux1299/.conda/envs/bad_mutations

# User defined input arguments
# Full path to a list of lists (to utilize GNU parallel and job arrays)
FASTA_LIST_OF_LISTS=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/align_lists/all_cds_hvulgare_list_of_lists.txt
# Full path to the config file
CONFIG_FILE=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/config.txt
# Full path to the out directory
# Note: subdirectories for each list of lists will automatically be generated (e.g., ${OUT_DIR}/list-00)
OUT_DIR=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/MSA_output
# Full path to the BAD_Mutations.py script
BAD_MUT_SCRIPT=~/Software/BAD_Mutations/BAD_Mutations.py
# Full path to where we want to store the log files output from parallel
LOG_FILE_DIR=${OUT_DIR}/all_parallel_log_files

#------------------------------
# Prepare array for list of lists
FASTA_LISTS_ARR=($(cat ${FASTA_LIST_OF_LISTS}))
# Determine maximum array limit
MAX_ARRAY_LIMIT=$[${#FASTA_LISTS_ARR[@]} - 1]
echo "Maximum array limit is ${MAX_ARRAY_LIMIT}."

# Get current FASTA list we are processing
CURR_FASTA_LIST=${FASTA_LISTS_ARR[${SLURM_ARRAY_TASK_ID}]}
# Generate array for current FASTA list we are processing
CURR_FASTA_ARR=($(cat ${CURR_FASTA_LIST}))

# Make subdirectory for current FASTA list we are processing (for easier tracking)
CURR_LIST_PREFIX=$(basename ${CURR_FASTA_LIST} .txt)

# Check if out directories exist, if not make them
mkdir -p ${OUT_DIR} \
    ${OUT_DIR}/${CURR_LIST_PREFIX} \
    ${OUT_DIR}/${CURR_LIST_PREFIX}/all_log_files \
    ${LOG_FILE_DIR}

# For experimenting with Slurm job arrays
echo "This is array index ${SLURM_ARRAY_TASK_ID}. Processing the list: ${CURR_FASTA_LIST}."

function bad_mut_align() {
    local bad_mut_script="$1"
    local config_file="$2"
    local curr_fasta_file="$3"
    local out_dir="$4"
    local curr_list_prefix="$5"
    # Generate log file prefix
    log_prefix=$(basename ${curr_fasta_file} .fa)
    set -x
    echo "Start processing: ${curr_fasta_file}"
    ${bad_mut_script} -v DEBUG align \
        -c ${config_file} \
        -f ${curr_fasta_file} \
        -o ${out_dir}/${curr_list_prefix} &> ${out_dir}/${curr_list_prefix}/all_log_files/${log_prefix}.log
    echo "Done: ${curr_fasta_file}"
    set +x
}

export -f bad_mut_align

# For each fasta list, parallelize across each fasta sequence record (one record per file in this list)
# Keep job log for parallel processes so upon resubmitting job, parallel can just re-run
#   samples that don't have an exit status of 0.
parallel --resume-failed --joblog ${LOG_FILE_DIR}/bad_mut_align.sh.${SLURM_ARRAY_TASK_ID}.log bad_mut_align ${BAD_MUT_SCRIPT} ${CONFIG_FILE} {} ${OUT_DIR} ${CURR_LIST_PREFIX} ::: ${CURR_FASTA_ARR[@]}
