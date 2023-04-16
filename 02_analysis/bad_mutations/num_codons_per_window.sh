#!/bin/bash

set -e
set -o pipefail

# Dependencies
module load bedtools/2.29.2

# User provided input arguments
# Indexed fasta file, .fai
GENOME_FILE="$1"
GFF="$2"
OUT_DIR="$3"
# Window size in bp
WIN_SIZE="$4"
# Feature we are pulling, must match pattern in column 6 in GFF file
FEATURE_TYPE="$5"

#------------
# Convert bp to Mbp
win_size_mbp=$(echo "${WIN_SIZE} / 1000000" | bc)
# Prepare output file column names
printf '#chr\tstart_pos\tend_pos\tfeature_type\tsum_length\tcount\tnum_codons\n' > ${OUT_DIR}/${FEATURE_TYPE}_${win_size_mbp}Mbp_num_codons_per_window.txt
# Prepare windows, pull feature type, sum and count, then get number of codons
bedtools makewindows -g ${GENOME_FILE} -w ${WIN_SIZE} | bedtools intersect -a - -b ${GFF} -wo | awk -v feature_type="${FEATURE_TYPE}" '{if ($6==feature_type) {print $0}}' | bedtools groupby -g 1,2,3,6 -c 13 -o sum,count | awk '{FS="\t"; OFS="\t"; print $1,$2,$3,$4,$5,$6,$5/3}' >> ${OUT_DIR}/${FEATURE_TYPE}_${win_size_mbp}Mbp_num_codons_per_window.txt
