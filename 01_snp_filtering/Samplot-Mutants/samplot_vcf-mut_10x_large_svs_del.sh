#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=8gb
#SBATCH --tmp=6gb
#SBATCH -t 00:30:00
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
source activate /panfs/jay/groups/9/morrellp/liux1299/.conda/envs/samplot_env

# User provided input arguments
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/samplot-mut_10x_large_svs_dels"
# VCF file
# All of the following VCFs have excluded HOM REF sites
DEL_VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/mut_3_lines_large_svs.DEL.private.callable.noMorexDiffs.vcf.gz"
# BAM file
#BAM1="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M01-3-3/outs/M01-3-3_phased_possorted_bam.bam"
BAM2="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M20-2-2/outs/M20-2-2_phased_possorted_bam.bam"
#BAM3="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M29-2-2/outs/M29-2-2_phased_possorted_bam.bam"
#SAMP1="M01"
SAMP2="M20"
#SAMP3="M29"
# Important regions
NUM_THREADS="16"

#----------------
# Make output directories
mkdir -p ${OUT_DIR} ${OUT_DIR}/DEL

# 10x Genomics VCF
# Note: Sample IDs have to come after list of BAM files to get it working
echo "Generating images for DEL..."
cd ${OUT_DIR}/DEL
samplot vcf \
    --plot_all \
    --vcf ${DEL_VCF} \
    -d ${OUT_DIR}/DEL/ \
    -O png \
    -a \
    -b ${BAM2} \
    -n ${SAMP2} \
    --threads ${NUM_THREADS} > ${OUT_DIR}/DEL/samplot_commands.sh
