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

# Use samplot to generate images of variants for each VCF individually

# Dependencies
module load python3/3.8.3_anaconda2020.07_mamba
source activate /panfs/roc/groups/9/morrellp/liux1299/.conda/envs/samplot_env

# User provided input arguments
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration/png_samplot_morex_vcfs"
# BAM files, one per variable
BAM1="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/morex-sample2/morex-sample2_phased_possorted_bam.bam"
BAM2="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_nanopore_V2_partsRef/Morex_nanopore_V2_partsRef_90_wRG.bam"
BAM3="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/pacbio_morex/pacbio_morex_v2/Morex_pacbio_90_wRG.bam"
# VCF files, one per variable in the same order as the BAM files
VCF1="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/morex-sample2/Filtered/morex-sample2_filtered_pass1.vcf"
VCF2="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_nanopore_V2_partsRef/Morex_nanopore_V2_partsRef_90.vcf"
VCF3="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/pacbio_morex/pacbio_morex_v2/Morex_pacbio_90.vcf"

#----------------
# Make output directories
mkdir -p ${OUT_DIR} ${OUT_DIR}/10xGenomics_vcf ${OUT_DIR}/ont_vcf ${OUT_DIR}/pacbio_vcf

# 10x Genomics VCF
echo "Generating images for 10x Genomics vcf..."
cd ${OUT_DIR}/10xGenomics_vcf
samplot vcf \
    --filter "SVTYPE == 'DEL'" \
    --filter "SVTYPE == 'INV'" \
    --filter "SVTYPE == 'BND'" \
    --filter "SVTYPE == 'DUP'" \
    --vcf ${VCF1} \
    -d ${OUT_DIR}/10xGenomics_vcf/ \
    -O png\
    -b ${BAM1} > samplot_commands.sh

# Nanopore VCF
echo "Generating images for Nanopore vcf..."
cd ${OUT_DIR}/ont_vcf
samplot vcf \
    --filter "SVTYPE == 'DEL'" \
    --filter "SVTYPE == 'INS'" \
    --filter "SVTYPE == 'INV'" \
    --filter "SVTYPE == 'BND'" \
    --filter "SVTYPE == 'DUP'" \
    --vcf ${VCF2} \
    -d ${OUT_DIR}/ont_vcf/ \
    -O png \
    -b ${BAM2} > samplot_commands.sh

# PacBio VCF
echo "Generating images for PacBio vcf..."
cd ${OUT_DIR}/pacbio_vcf
samplot vcf \
    --filter "SVTYPE == 'DEL'" \
    --filter "SVTYPE == 'INS'" \
    --filter "SVTYPE == 'INV'" \
    --filter "SVTYPE == 'BND'" \
    --filter "SVTYPE == 'DUP'" \
    --vcf ${VCF3} \
    -d ${OUT_DIR}/pacbio_vcf/ \
    -O png \
    -b ${BAM3} > samplot_commands.sh
