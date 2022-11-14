#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=12gb
#SBATCH --tmp=10gb
#SBATCH -t 00:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

set -e
set -o pipefail

# Dependencies
# Mosdepth version 0.3.1
export PATH=${PATH}:/panfs/jay/groups/9/morrellp/shared/Software/mosdepth

# Input arguments
BAM="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/M01_ont_partsRefv3/M01_ont_partsRefv3_90_wRG.bam"
SAMPLE_NAME="M01_ont"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/mosdepth_coverage"

#-----------------------
# Go into output directory
cd ${OUT_DIR}

# by setting these ENV vars, we can control the output labels (4th column)
export MOSDEPTH_Q0=NO_COVERAGE   # 0 -- defined by the arguments to --quantize
export MOSDEPTH_Q1=LOW_COVERAGE  # 1..4
export MOSDEPTH_Q2=CALLABLE      # 5..149
export MOSDEPTH_Q3=HIGH_COVERAGE # 150 ...

# Get BED with labels for coverage
mosdepth -n --quantize 0:1:5:150: ${SAMPLE_NAME} ${BAM}
