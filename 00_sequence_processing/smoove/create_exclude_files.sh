#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=6gb
#SBATCH --tmp=4gb
#SBATCH -t 00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

set -e
set -o pipefail

# Create sample specific exclude files utilizing Slurm job arrays

# Dependencies
module load bedtools/2.29.2

# User provided input arguments
# List of high coverage bed files
BED_LIST="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/smoove_processing/masked_regions/per_sample_high_depth/high_cov_bed_list.txt"
# Regions where there are stretches of N's in the reference genome
GAP_BED="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/stretches_of_Ns/Barley_MorexV3_pseudomolecules_parts_missing.bed"
# Where to store output files
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/smoove_processing/masked_regions/per_sample_gap_and_high_depth"

#----------------
mkdir -p ${OUT_DIR}

# Prepare array for utilizing Slurm job arrays
BED_ARR=($(cat ${BED_LIST}))
# Determine maximum array limit
MAX_ARRAY_LIMIT=$[${#BED_ARR[@]} - 1]
echo "Maximum array limit is ${MAX_ARRAY_LIMIT}."

# Get the current BED file we are processing
CURR_BED=${BED_ARR[${SLURM_ARRAY_TASK_ID}]}
echo "Currently processing BED file: ${CURR_BED}"

function concat_and_merge_bed() {
    local gap_bed="$1"
    local high_cov_bed="$2"
    local out_dir="$3"
    # Get sample name
    sample_name=$(basename ${high_cov_bed} .high_cov.bed)
    # Concat gap BED and high cov BED, then merge near regions
    cat ${gap_bed} <(cut -f1-3 ${high_cov_bed}) | bedtools sort -i stdin | bedtools merge -d 10 -i stdin > ${out_dir}/${sample_name}.exclude.bed
}

export -f concat_and_merge_bed

concat_and_merge_bed ${GAP_BED} ${CURR_BED} ${OUT_DIR}
