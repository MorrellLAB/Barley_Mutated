#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=8gb
#SBATCH --tmp=6gb
#SBATCH -t 05:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Add @RG header to BAM file

# Dependencies
module load samtools/1.9

# User provided input arguments
BAM="/scratch.global/liux1299/nanopore_morex_v3/Morex_ont_partsRefv3/Morex_ont_partsRefv3_90.bam"
RG_LINE="@RG\tID:Morex_ont\tPL:ONT\tSM:Morex_ont"
OUT_DIR="/scratch.global/liux1299/nanopore_morex_v3/Morex_ont_partsRefv3"
NUM_THREADS="16"

#-----------------
mkdir -p ${OUT_DIR}

OUT_PREFIX=$(basename ${BAM} .bam)

# Add @RG header line to bam
samtools addreplacerg -r ${RG_LINE} --threads ${NUM_THREADS} -o ${OUT_DIR}/${OUT_PREFIX}_wRG.bam ${BAM}
# Index SAM
samtools index -b ${OUT_DIR}/${OUT_PREFIX}_wRG.bam
