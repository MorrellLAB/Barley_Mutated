#!/bin/bash

# Generate VCF of positions that overlap with low complexity regions
# For exploration and filter tuning purposes, output will be file with "temp_" in filename
# These positions can be used to jump to locations when exploring SVs in IGV

# Dependencies
module load bedtools/2.29.2

# User provided input arguments
# Command line input
# Input VCF
VCF=$(realpath "$1")
OUT_DIR="$2"

# Predefined arguments, can add more if needed
# Low complexity BED file(s) with varying entropy levels
LOW_COMP1="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/entropy_masked/Barley_MorexV3_pseudomolecules_parts.entropy_0.9_masked.bed"
LOW_COMP2="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/entropy_masked/Barley_MorexV3_pseudomolecules_parts.entropy_0.8_masked.bed"
LOW_COMP3="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/entropy_masked/Barley_MorexV3_pseudomolecules_parts.entropy_0.7_masked.bed"
# Associated entropy values that match order defined in LOW_COMP1, LOW_COMP2, etc.
# These will be used in the output filenames
ENTROPY1="0.9"
ENTROPY2="0.8"
ENTROPY3="0.7"

#------------------------
function check_filter_stringency() {
    local vcf_file="$1"
    local num_sites="$2"
    if [[ ${num_sites} == 0 ]]; then
        echo "No sites left after filtering. Try using less stringent criteria. File with no sites is: ${vcf_file}. Exiting..." >&2
        exit 8 # If not sites left, error out with message
    fi
}

export -f check_filter_stringency

function count_sites() {
    local vcf="$1"
    local log_file="$2"
    # Get the number of sites left after filtering
    if [[ "${vcf}" == *"gz"* ]]; then
        # We are working with gzipped VCF file
        num_sites=$(zgrep -v "#" ${vcf} | wc -l)
    else
        # We are working with uncompressed vcf
        num_sites=$(grep -v "#" ${vcf} | wc -l)
    fi
    # Append the number of sites remaining to file
    printf "${vcf}\t${num_sites}\n" >> ${log_file}
    check_filter_stringency ${vcf} ${num_sites}
}

export -f count_sites

function pull_low_complex_overlaps() {
    local vcf="$1"
    local low_complex_bed="$2"
    local entropy_val="$3"
    local out_dir="$4"
    # Prepare header
    if [[ "${vcf}" == *"vcf.gz" ]]; then
        # Get vcf basename
        vcf_bn=$(basename ${vcf} .vcf.gz)
        # Working with gzipped vcf file
        zgrep "#" ${vcf} > ${out_dir}/temp_entropy_${entropy_val}_x_${vcf_bn}.vcf
    elif [[ "${vcf}" == *"vcf" ]]; then
        # Get vcf basename
        vcf_bn=$(basename ${vcf} .vcf)
        # Working with uncompressed vcf
        grep "#" ${vcf} > ${out_dir}/temp_entropy_${entropy_val}_x_${vcf_bn}.vcf
    fi
    # Intersect VCF with low complexity BED file
    bedtools intersect -wa -a ${vcf} -b ${low_complex_bed} | sort -uV -k1,1 -k2,2 >> ${out_dir}/temp_entropy_${entropy_val}_x_${vcf_bn}.vcf
    # Count the number of sites
    count_sites ${out_dir}/temp_entropy_${entropy_val}_x_${vcf_bn}.vcf ${OUT_DIR}/low_complexity_num_sites.log
}

export -f pull_low_complex_overlaps

# Make output directory if it doesn't exist
mkdir -p ${OUT_DIR}

# Get starting number of sites
count_sites ${VCF} ${OUT_DIR}/low_complexity_num_sites.log

# Pull SVs that overlap low complexity region for each entropy value
echo "Running for entropy value: ${ENTROPY1}..."
pull_low_complex_overlaps ${VCF} ${LOW_COMP1} ${ENTROPY1} ${OUT_DIR}

echo "Running for entropy value: ${ENTROPY2}..."
pull_low_complex_overlaps ${VCF} ${LOW_COMP2} ${ENTROPY2} ${OUT_DIR}

echo "Running for entropy value: ${ENTROPY3}..."
pull_low_complex_overlaps ${VCF} ${LOW_COMP3} ${ENTROPY3} ${OUT_DIR}
