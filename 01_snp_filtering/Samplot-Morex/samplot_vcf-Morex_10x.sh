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
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration_morex_v3/png_samplot_morex_vcfs/morex_10x"
# VCF file
# All of the following VCFs have excluded HOM REF sites
DEL_VCF="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/morex-sample2_dels.filtered_pass1.SVs_not_BND.vcf"
LSV_VCF="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/morex-sample2_large_svs.filtered_pass1.SVs_not_BND.vcf"
PHV_VCF="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/morex-sample2_phased_variants.dp_ann_filt.vcf"
# BAM file
BAM="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/outs/morex-sample2_phased_possorted_bam.bam"
# Important regions
NUM_THREADS="16"

#----------------
# Make output directories
mkdir -p ${OUT_DIR} ${OUT_DIR}/DEL ${OUT_DIR}/Large_SVs ${OUT_DIR}/Phased_Variants

# 10x Genomics VCF
echo "Generating images for DEL..."
cd ${OUT_DIR}/DEL
samplot vcf \
    --plot_all \
    --vcf ${DEL_VCF} \
    -d ${OUT_DIR}/DEL/ \
    -O png \
    -b ${BAM} \
    --threads ${NUM_THREADS} > ${OUT_DIR}/DEL/samplot_commands.sh

echo "Generating images for Large_SVs..."
cd ${OUT_DIR}/Large_SVs
samplot vcf \
    --plot_all \
    --vcf ${LSV_VCF} \
    -d ${OUT_DIR}/Large_SVs/ \
    -O png \
    -b ${BAM} \
    --threads ${NUM_THREADS} > ${OUT_DIR}/Large_SVs/samplot_commands.sh

echo "Generating images for Phased Variants..."
cd ${OUT_DIR}/Phased_Variants
samplot vcf \
    --plot_all \
    --vcf ${PHV_VCF} \
    -d ${OUT_DIR}/Phased_Variants/ \
    -O png \
    -b ${BAM} \
    --threads ${NUM_THREADS} > ${OUT_DIR}/Phased_Variants/samplot_commands.sh
