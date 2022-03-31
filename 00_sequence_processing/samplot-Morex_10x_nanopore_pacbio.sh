#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=36gb
#SBATCH --tmp=22gb
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Use samplot to generate images of variants using GNU parallel

# Dependencies
module load parallel/20210822
module load python3/3.8.3_anaconda2020.07_mamba
source activate /panfs/roc/groups/9/morrellp/liux1299/.conda/envs/samplot_env

# User provided input arguments
REGIONS="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration/morex_10xGenomics_nanopore_pacbio_merged.bed"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration/png_samplot_morex_10xGenomics_nanopore_pacbio"
# Sample IDs, one per variable
# Must be in the same order as the BAM files
SAMPLE_ID1="Morex_10x"
SAMPLE_ID2="Morex_nanopore"
SAMPLE_ID3="Morex_pacbio"
# BAM files, one per variable
BAM1="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/morex-sample2/morex-sample2_phased_possorted_bam.bam"
BAM2="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_nanopore_V2_partsRef/Morex_nanopore_V2_partsRef_90_wRG.bam"
BAM3="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/pacbio_morex/pacbio_morex_v2/Morex_pacbio_90_wRG.bam"
# This bed file lists all intervals with 'Ns' in the reference genome
ANNOTATE="/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Barley_Morex_V2_pseudomolecules_parts_missing.bed.gz"
# Type of variant
VAR_TYPE="DEL"

#----------------
mkdir -p ${OUT_DIR}

function run_samplot() {
    local regions="$1"
    local out_dir="$2"
    local annotations="$3"
    local var_type="$4"
    local curr_region_line="$5"
    local samp_id1="$6"
    local samp_id2="$7"
    local samp_id3="$8"
    local bam1="$9"
    local bam2="${10}"
    local bam3="${11}"
    curr_chr=$(echo "${curr_region_line}" | cut -d ';' -f 1)
    curr_start=$(echo "${curr_region_line}" | cut -d ';' -f 2)
    curr_end=$(echo "${curr_region_line}" | cut -d ';' -f 3)
    # Run samplot and generate image for region
    samplot plot \
        --important_regions "${regions}" \
        --sample_ids ${samp_id1} ${samp_id2} ${samp_id3} \
        --out-dir "${out_dir}/" \
        -O png \
        -b ${bam1} ${bam2} ${bam3} \
        -A "${annotations}" \
        -t "${var_type}" \
        --chrom "${curr_chr}" \
        --start "${curr_start}" \
        --end "${curr_end}" \
        --output_file "${out_dir}/${var_type}_${curr_chr}_${curr_start}_${curr_end}.png"
}

export -f run_samplot

# Reformat regions line and store in array
bed_arr=($(sed -e 's,\t,;,g' ${REGIONS}))

parallel --verbose run_samplot "${REGIONS}" "${OUT_DIR}" "${ANNOTATE}" "${VAR_TYPE}" {} ${SAMPLE_ID1} ${SAMPLE_ID2} ${SAMPLE_ID3} ${BAM1} ${BAM2} ${BAM3} ::: "${bed_arr[@]}"
