#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=12gb
#SBATCH --tmp=10gb
#SBATCH -t 01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

set -e
set -o pipefail

# Dependencies
module load bedtools/2.29.2
# Mosdepth version 0.3.1
export PATH=${PATH}:/panfs/jay/groups/9/morrellp/shared/Software/mosdepth

# Input arguments
BAM="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/outs/morex-sample2_phased_possorted_bam.bam"
SAMPLE_NAME="morex-sample2"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/mosdepth_coverage"
WIN_SIZE="10"
LOW_COV_CUTOFF="5"
NUM_THREADS="4"

#-----------------------
# Make output directory
mkdir -p ${OUT_DIR}
# Go into output directory
cd ${OUT_DIR}

# by setting these ENV vars, we can control the output labels (4th column)
export MOSDEPTH_Q0=NO_COVERAGE   # 0 -- defined by the arguments to --quantize
export MOSDEPTH_Q1=LOW_COVERAGE  # 1..4
export MOSDEPTH_Q2=CALLABLE      # 5..149
export MOSDEPTH_Q3=HIGH_COVERAGE # 150 ...

# Get BED with labels for coverage
mosdepth --by ${WIN_SIZE} --threads ${NUM_THREADS} --fast-mode --no-per-base --quantize 0:1:5:150: ${SAMPLE_NAME} ${BAM}

# Use mosdepth output file
# Remove chrUn and merge overlapping or "book-ended' features
zcat ${OUT_DIR}/${SAMPLE_NAME}.regions.bed.gz | awk -v my_cutoff=${LOW_COV_CUTOFF} '$4 < my_cutoff' | grep -v "chrUn" | bedtools merge -i - > ${OUT_DIR}/${SAMPLE_NAME}.regions.lt${LOW_COV_CUTOFF}.bed
