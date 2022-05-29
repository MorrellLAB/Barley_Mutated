#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=20gb
#SBATCH --tmp=16gb
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Use samplot to generate images of variants for each VCF individually
# Works with single VCF file!

# Dependencies
module load python3/3.8.3_anaconda2020.07_mamba
source activate /panfs/roc/groups/9/morrellp/liux1299/.conda/envs/samplot_env

# User provided input arguments
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration_morex_v3/png_samplot_morex_vcfs/morex_10x_ont_pacbio"
# VCF file
VCF="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration_morex_v3/morex_10x_ont_pacbio_noHomRef.vcf.gz"
# BAM files, one per variable
# In the same order as the samples in the VCF
BAM_10x="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/outs/morex-sample2_phased_possorted_bam.bam"
BAM_ONT="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3/Morex_ont_partsRefv3_90_wRG.bam"
BAM_PACBIO="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/pacbio_morex/pacbio_morex_v3/Morex_pacbio_90_wRG.bam"
NUM_THREADS="16"

#----------------
# Make output directories
mkdir -p ${OUT_DIR} ${OUT_DIR}/DEL ${OUT_DIR}/INS ${OUT_DIR}/INV ${OUT_DIR}/DUP

# DEL
echo "Generating images for DEL..."
cd ${OUT_DIR}/DEL
samplot vcf \
    --filter "SVTYPE == 'DEL'" \
    --vcf ${VCF} \
    -d ${OUT_DIR}/DEL/ \
    -O png \
    -b ${BAM_10x} ${BAM_ONT} ${BAM_PACBIO} \
    --threads ${NUM_THREADS} > ${OUT_DIR}/DEL/samplot_commands.sh

# INS
echo "Generating images for INS..."
cd ${OUT_DIR}/INS
samplot vcf \
    --filter "SVTYPE == 'INS'" \
    --vcf ${VCF} \
    -d ${OUT_DIR}/INS/ \
    -O png \
    -b ${BAM_10x} ${BAM_ONT} ${BAM_PACBIO} \
    --threads ${NUM_THREADS} > ${OUT_DIR}/INS/samplot_commands.sh

# INV
echo "Generating images for INV..."
cd ${OUT_DIR}/INV
samplot vcf \
    --filter "SVTYPE == 'INV'" \
    --vcf ${VCF} \
    -d ${OUT_DIR}/INV/ \
    -O png \
    -b ${BAM_10x} ${BAM_ONT} ${BAM_PACBIO} \
    --threads ${NUM_THREADS} > ${OUT_DIR}/INV/samplot_commands.sh

# DUP
echo "Generating images for DUP..."
cd ${OUT_DIR}/DUP
samplot vcf \
    --filter "SVTYPE == 'DUP'" \
    --vcf ${VCF} \
    -d ${OUT_DIR}/DUP/ \
    -O png \
    -b ${BAM_10x} ${BAM_ONT} ${BAM_PACBIO} \
    --threads ${NUM_THREADS} > ${OUT_DIR}/DUP/samplot_commands.sh
