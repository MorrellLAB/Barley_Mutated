#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=20gb
#SBATCH --tmp=16gb
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

set -e
set -o pipefail

# Use samplot to generate images of variants for each VCF individually
# Works with single VCF file!
# Utilizes Slurm job arrays, one per variant type

# Dependencies
module load python3/3.8.3_anaconda2020.07_mamba
source activate /panfs/roc/groups/9/morrellp/liux1299/.conda/envs/samplot_env

# User provided input arguments
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration_morex_v3/png_samplot_morex_vcfs/smoove_morex10x_wgs_mut_hybrid"
# VCF file
VCF="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/smoove_processing/results_smoove/results_genotyped/mut_barley_cohort_sites.smoove.square.anno.vcf.gz"
# BAM list, one per line sorted in the same order as the samples in the VCF
# Here's how the sorting was done:
#   In dir: ~/Projects/Mutant_Barley/smoove_processing/results_smoove/results_genotyped
#   cat /panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/smoove_processing/bam_symlinks/wgs_bam_list.txt | sed 's,\(.*\)/,\1\t,' | sort -k2,2 | sed 's,\t,/,' > wgs_bam_list_vcf_sample_order.txt
BAM_LIST="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/smoove_processing/results_smoove/results_genotyped/wgs_bam_list_vcf_sample_order.txt"
# BAM files, one per variable sorted in the
# In the same order as the samples in the VCF
#BAM_10x="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/outs/morex-sample2_phased_possorted_bam.bam"
#BAM_ONT="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3/Morex_ont_partsRefv3_90_wRG.bam"
#BAM_PACBIO="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/pacbio_morex/pacbio_morex_v3/Morex_pacbio_90_wRG.bam"
NUM_THREADS="16"

# Type of variants to plot
VAR_TYPE1="DEL"
VAR_TYPE2="INS"
VAR_TYPE3="INV"
VAR_TYPE4="DUP"
VAR_TYPE5="BND"

#----------------
# Make output directories
mkdir -p ${OUT_DIR} ${OUT_DIR}/${VAR_TYPE1} ${OUT_DIR}/${VAR_TYPE2} ${OUT_DIR}/${VAR_TYPE3} ${OUT_DIR}/${VAR_TYPE4} ${OUT_DIR}/${VAR_TYPE5}

# Store variant types in bash array
VAR_ARR=("${VAR_TYPE1}" "${VAR_TYPE2}" "${VAR_TYPE3}" "${VAR_TYPE4}" "${VAR_TYPE5}")

# Prepare to run as Slurm job array
# Determine maximum array limit
MAX_ARRAY_LIMIT=$[${#VAR_ARR[@]} - 1]
echo "Maximum array limit is ${MAX_ARRAY_LIMIT}."
# Get current variant type for job array
CURR_VAR_TYPE=${VAR_ARR[${SLURM_ARRAY_TASK_ID}]}
echo "Current variant type is: ${CURR_VAR_TYPE}"

# Run samplot vcf for current variant type
echo "Generating images for ${CURR_VAR_TYPE}..."
cd ${OUT_DIR}/${CURR_VAR_TYPE}
samplot vcf \
    --filter "SVTYPE == ${CURR_VAR_TYPE}" \
    --vcf ${VCF} \
    -d ${OUT_DIR}/${CURR_VAR_TYPE}/ \
    -O png \
    -b $(cat ${BAM_LIST} | tr '\n' ' ') \
    --threads ${NUM_THREADS} > ${OUT_DIR}/${CURR_VAR_TYPE}/samplot_commands.sh
