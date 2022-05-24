#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=48gb
#SBATCH --tmp=36gb
#SBATCH -t 01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# This script calls on the extract_high_depth_regions.py script to
#   extract high coverage regions that could be problematic when
#   running smoove.

# Dependencies
module load parallel/20210822
module load python3/3.8.3_anaconda2020.07_mamba
# Python script: extract_high_depth_regions.py
export PATH=${PATH}:/panfs/roc/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/00_sequence_processing/smoove

# User provided input arguments
# List of filepaths to per sample regions.bed.gz files
BED_LIST="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/coverage_wgs/mosdepth_bed_list.txt"
# For barley, we'll exclude chrUn when generating the high coverage bed list
EXCLUDE_CHR="chrUn"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/smoove_processing/masked_regions"

#--------------------
# Make output directories
mkdir -p ${OUT_DIR} ${OUT_DIR}/per_sample_high_depth

function extract_high_depth_regions() {
    local bedgz_file="$1"
    local exclude_chr="$2"
    local out_dir="$3"
    # Get sample name of current bed file
    sample_name=$(basename ${bedgz_file} .regions.bed.gz)
    # Get directory of current bed file
    bed_dir=$(dirname ${bedgz_file})
    # Exclude chromosome specified
    zgrep -v "${exclude_chr}" ${bedgz_file} | gzip -c > ${bed_dir}/${sample_name}_no${exclude_chr}.regions.bed.gz
    # Extract high depth regions
    extract_high_depth_regions.py ${bed_dir}/${sample_name}_no${exclude_chr}.regions.bed.gz ${sample_name} ${out_dir}
}

export -f extract_high_depth_regions

parallel --verbose extract_high_depth_regions {} ${EXCLUDE_CHR} ${OUT_DIR}/per_sample_high_depth :::: ${BED_LIST}
