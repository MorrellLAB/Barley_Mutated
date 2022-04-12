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

# Use samplot to generate images of variants using GNU parallel and Slurm job arrays
# Usage: sbatch --array=0-1 samplot-Morex_ont_pacbio.sh
# Where: "0-1" corresponds to the number of variant types stored in a bash array

# Dependencies
module load parallel/20210822
module load python3/3.8.3_anaconda2020.07_mamba
source activate /panfs/roc/groups/9/morrellp/liux1299/.conda/envs/samplot_env

# User provided input arguments
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration_morex_v3/png_samplot_morex_vcfs/ont_pacbio_plot"
# Sample IDs, one per variable
# Must be in the same order as the BAM files
SAMPLE_ID1="morex_ont"
SAMPLE_ID2="morex_pacbio"
# BAM files, one per variable
BAM1="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3/Morex_ont_partsRefv3_90_wRG.bam"
BAM2="/scratch.global/liux1299/sra_pacbio/pacbio_morex_v3/Morex_pacbio/Morex_pacbio_90_wRG.bam"
# VCF files, one per variable in the same order as the BAM files
#VCF1="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3/filtered/morex_ont_noHomRef.vcf"
#VCF2="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/pacbio_morex/pacbio_morex_v3/filtered/morex_pacbio_noHomRef.vcf"
# This bed file lists all intervals with 'Ns' in the reference genome
ANNOTATE1="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/stretches_of_Ns/Barley_MorexV3_pseudomolecules_parts_missing.bed.gz"
# This bed file lists all intervals with TEs
ANNOTATE2="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/repeat_annotation/TEanno-v1.0__200416_MorexV3_pseudomolecules.sorted.parts.minimal.bed.gz"
# Type of variant
VAR_TYPE1="DEL"
VAR_TYPE2="INS"
VAR_TYPE3="INV"
#VAR_TYPE4="BND"
# Regions file for each variant type, in the same order as type of variant
REGION1="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration_morex_v3/DEL_morex_ONT_PacBio_sorted.bed"
REGION2="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration_morex_v3/INS_morex_ONT_PacBio_sorted.bed"
REGION3="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration_morex_v3/INV_morex_ONT_PacBio_sorted.bed"

#----------------
# Store variant types in bash array
#VAR_ARR=("${VAR_TYPE1}" "${VAR_TYPE2}" "${VAR_TYPE3}" "${VAR_TYPE4}")
VAR_ARR=("${VAR_TYPE1}" "${VAR_TYPE2}" "${VAR_TYPE3}")
REGION_ARR=("${REGION1}" "${REGION2}" "${REGION3}")
# Make output directories
mkdir -p ${OUT_DIR}
# Make output subdirectories
for i in ${VAR_ARR[@]}
do
    mkdir -p ${OUT_DIR}/${i}
done

# Prepare to run as Slurm job array
# Determine maximum array limit
MAX_ARRAY_LIMIT=$[${#VAR_ARR[@]} - 1]
echo "Maximum array limit is ${MAX_ARRAY_LIMIT}."
# Get current variant type for job array
CURR_VAR_TYPE=${VAR_ARR[${SLURM_ARRAY_TASK_ID}]}
# Get current regions file based on the variant type for job array
REGIONS=${REGION_ARR[${SLURM_ARRAY_TASK_ID}]}
echo "Current variant type is: ${CURR_VAR_TYPE}"
echo "Current regions file is: ${REGIONS}"

function run_samplot() {
    local regions="$1"
    local out_dir="$2"
    local annotations1="$3"
    local annotations2="$4"
    local var_type="$5"
    local curr_region_line="$6"
    local samp_id1="$7"
    local samp_id2="$8"
    local bam1="$9"
    local bam2="${10}"
    curr_chr=$(echo "${curr_region_line}" | cut -d ';' -f 1)
    curr_start=$(echo "${curr_region_line}" | cut -d ';' -f 2)
    curr_end=$(echo "${curr_region_line}" | cut -d ';' -f 3)
    # Run samplot and generate image for region
    samplot plot \
        --important_regions "${regions}" \
        --sample_ids ${samp_id1} ${samp_id2} \
        --out-dir "${out_dir}/" \
        -O png \
        -b ${bam1} ${bam2} \
        -A "${annotations1}" "${annotations2}" \
        -t "${var_type}" \
        --chrom "${curr_chr}" \
        --start "${curr_start}" \
        --end "${curr_end}" \
        --output_file "${out_dir}/${var_type}_${curr_chr}_${curr_start}_${curr_end}.png"
}

export -f run_samplot

# Reformat regions line and store in array
bed_arr=($(sed -e 's,\t,;,g' ${REGIONS}))
# Get the number of regions
total_lines=$(echo ${#bed_arr[@]})
# If total number of lines exceeds 53,083, we'll need to split the
#   bed_arr into batches otherwise we'll get the following error from GNU parallel:
#   bash: /panfs/roc/msisoft/parallel/20210822/bin/parallel: Argument list too long
# We'll use 50000 for a nicer number
if [ "${total_lines}" -gt 50000 ]; then
    echo "Total number of regions exceeds 50000. Splitting GNU parallel job into batches of regions."
    # Split bed_arr into batches
    bed_arr_batch1=($(sed -n "1,50000p" ${REGIONS} | sed -e 's,\t,;,g'))
    bed_arr_batch2=($(sed -n "50001,${total_lines}p" ${REGIONS} | sed -e 's,\t,;,g'))
    # Check the number of lines in each batch
    echo "Batch1: ${#bed_arr_batch1[@]}"
    echo "Batch2: ${#bed_arr_batch2[@]}"
    # Current variant type, batch1
    parallel --verbose run_samplot "${REGIONS}" "${OUT_DIR}/${CURR_VAR_TYPE}" "${ANNOTATE1}" "${ANNOTATE2}" "${CURR_VAR_TYPE}" {} ${SAMPLE_ID1} ${SAMPLE_ID2} ${BAM1} ${BAM2} ::: "${bed_arr_batch1[@]}"
    # Current variant type, batch2
    parallel --verbose run_samplot "${REGIONS}" "${OUT_DIR}/${CURR_VAR_TYPE}" "${ANNOTATE}" "${CURR_VAR_TYPE}" {} ${SAMPLE_ID1} ${SAMPLE_ID2} ${BAM1} ${BAM2} ::: "${bed_arr_batch2[@]}"
else
    echo "Total number of regions doesn't exceed 50000. Running GNU parallel job as a single batch."
    # Run as a single batch
    # Current variant type
    parallel --verbose run_samplot "${REGIONS}" "${OUT_DIR}/${CURR_VAR_TYPE}" "${ANNOTATE1}" "${ANNOTATE2}" "${CURR_VAR_TYPE}" {} ${SAMPLE_ID1} ${SAMPLE_ID2} ${BAM1} ${BAM2} ::: "${bed_arr[@]}"
fi
