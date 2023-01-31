#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=8gb
#SBATCH --tmp=6gb
#SBATCH -t 03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Run cuteSV on ONT data

# Dependencies
module load python3/3.8.3_anaconda2020.07_mamba
# Conda env
source activate /panfs/jay/groups/9/morrellp/liux1299/.conda/envs/my_cutesv_env

# User provided input arguments
BAM="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/M29_ont_partsRefv3/M29_ont_partsRefv3_90_wRG.bam"
REF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta"
SAMPLE="M29_ont"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/cutesv_calls-${SAMPLE}"

#----------------
mkdir -p ${OUT_DIR}

cd ${OUT_DIR}
# ONT data settings
cuteSV ${BAM} \
    ${REF} \
    ${SAMPLE}.vcf \
    ${OUT_DIR} \
    --max_cluster_bias_INS 100 \
    --diff_ratio_merging_INS 0.3 \
    --max_cluster_bias_DEL 100 \
    --diff_ratio_merging_DEL 0.3 \
    --genotype
