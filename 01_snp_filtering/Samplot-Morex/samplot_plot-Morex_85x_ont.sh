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
module load parallel/20210822

# User provided input arguments
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/samplot-morex_85x_ont"
# Important regions
REGIONS_DEL="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_85x_ont/filtered/pseudo_pos/morex_85x_ont.noHomRef.geSup5.callable.DEL.bed.gz"
REGIONS_INS="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_85x_ont/filtered/pseudo_pos/morex_85x_ont.noHomRef.geSup5.callable.INS.bed.gz"
# BAM file (requires @RG header lines)
BAM="/scratch.global/liux1299/temp_ipk_morex_85x_ont/Morex_85X_sorted_wRG.bam"
SAMPLE_NAME="morex_85x_ont"

# BED file for annotations
# Uncallable regions
UNCALLABLE="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/uncallable_regions/morex_v3_combined_uncallable.low_complexity.nochrUn.pseudo_pos.bed.gz"

#----------------
# Make output directories
mkdir -p ${OUT_DIR} ${OUT_DIR}/DEL ${OUT_DIR}/INS

function run_samplot() {
    local regions="$1"
    local out_dir="$2"
    local annotations1="$3"
    local var_type="$4"
    local curr_region_line="$5"
    local samp_id1="$6"
    local bam1="$7"
    curr_chr=$(echo "${curr_region_line}" | cut -d ';' -f 1)
    curr_start=$(echo "${curr_region_line}" | cut -d ';' -f 2)
    curr_end=$(echo "${curr_region_line}" | cut -d ';' -f 3)
    # Run samplot and generate image for region
    # Note: -a flag is required for running SV-plaudit downstream
    samplot plot \
        --important_regions "${regions}" \
        --sample_ids ${samp_id1} \
        --out-dir "${out_dir}/" \
        -O png \
        -b ${bam1} \
        -A "${annotations1}" \
        -t "${var_type}" \
        --chrom "${curr_chr}" \
        --start "${curr_start}" \
        --end "${curr_end}" \
        -a \
        --output_file "${out_dir}/${var_type}_${curr_chr}_${curr_start}_${curr_end}.png"
}

export -f run_samplot

# Reformat regions line and store in array
if [[ "${REGIONS_DEL}" == *".gz" ]]; then
    del_bed_arr=($(zcat ${REGIONS_DEL} | sed -e 's,\t,;,g'))
else
    del_bed_arr=($(sed -e 's,\t,;,g' ${REGIONS_DEL}))
fi
# Get the number of regions
total_lines_del=$(echo ${#del_bed_arr[@]})

if [[ "${REGIONS_INS}" == *".gz" ]]; then
    ins_bed_arr=($(zcat ${REGIONS_INS} | sed -e 's,\t,;,g'))
else
    ins_bed_arr=($(sed -e 's,\t,;,g' ${REGIONS_INS}))
fi
# Get the number of regions
total_lines_ins=$(echo ${#ins_bed_arr[@]})

echo "Generating images for DEL..."
# Current variant type
parallel --verbose run_samplot "${REGIONS_DEL}" "${OUT_DIR}/DEL" "${UNCALLABLE}" "DEL" {} ${SAMPLE_NAME} ${BAM} ::: "${del_bed_arr[@]}"

echo "Generating images for INS..."
parallel --verbose run_samplot "${REGIONS_INS}" "${OUT_DIR}/INS" "${UNCALLABLE}" "INS" {} ${SAMPLE_NAME} ${BAM} ::: "${ins_bed_arr[@]}"
