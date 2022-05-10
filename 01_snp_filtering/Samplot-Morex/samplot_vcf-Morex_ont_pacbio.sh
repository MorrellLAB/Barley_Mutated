#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=20gb
#SBATCH --tmp=16gb
#SBATCH -t 36:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# DEPRECATED! DO NOT USE! For some reason, samplot doesn't work well with VCF files
# Use samplot to generate images of variants for each VCF individually

# Dependencies
module load python3/3.8.3_anaconda2020.07_mamba
source activate /panfs/roc/groups/9/morrellp/liux1299/.conda/envs/samplot_env

# User provided input arguments
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration_morex_v3/png_samplot_morex_vcfs/ont_pacbio_vcf"
# VCF file
VCF="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration_morex_v3/morex_ont_pacbio_noHomRef.vcf.gz"
# BAM files, one per variable
# In the same order as the samples in the VCF
BAM1="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3/Morex_ont_partsRefv3_90_wRG.bam"
BAM2="/scratch.global/liux1299/sra_pacbio/pacbio_morex_v3/Morex_pacbio/Morex_pacbio_90_wRG.bam"
# Important regions
REGIONS="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration_morex_v3/morex_ont_pacbio_noHomRef_noBND.bed"
NUM_THREADS="16"

#----------------
# Make output directories
mkdir -p ${OUT_DIR}

# 10x Genomics VCF
echo "Generating images for DEL..."
cd ${OUT_DIR}
samplot vcf \
    --plot_all \
    --vcf ${VCF} \
    -d ${OUT_DIR}/\
    -O png \
    --important_regions ${REGIONS} \
    -b ${BAM1} ${BAM2} \
    --threads ${NUM_THREADS} > ${OUT_DIR}/samplot_commands.sh
