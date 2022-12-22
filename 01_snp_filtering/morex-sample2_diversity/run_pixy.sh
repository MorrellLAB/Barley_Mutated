#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=12gb
#SBATCH --tmp=8gb
#SBATCH -t 06:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Dependencies
module load python3/3.8.3_anaconda2020.07_mamba
# Activate conda environment
source activate /panfs/jay/groups/9/morrellp/liux1299/.conda/envs/pixy_env

VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/pixy_input/morex-sample2.vcf.gz"
# Populations .txt file
POP_FILE="/panfs/jay/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/01_snp_filtering/morex-sample2_diversity/pop_file.txt"
# Window size (bp)
WIN_SIZE="100"
# Number of cores
NUM_CORES="16"
# Output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/pixy_output"

#---------------------
mkdir -p ${OUT_DIR}

cd ${OUT_DIR}

pixy --stats pi \
     --vcf ${VCF} \
     --populations ${POP_FILE} \
     --bypass_invariant_check 'yes' \
     --window_size ${WIN_SIZE} \
     --n_cores ${NUM_CORES}
