#!/bin/bash

set -e
set -o pipefail

# Dependencies
module load bedtools/2.29.2

# User provided input arguments
# Indexed fasta file, .fai
GENOME_FILE="$1"
# Extended BED file format with column 4 containing "Deleterious" vs "Tolerated" annotations
BED="$2"
OUT_DIR="$3"
OUT_PREFIX="$4"
# Window size in bp
WIN_SIZE="$5"

#------------
# Convert bp to Mbp
win_size_mbp=$(echo "${WIN_SIZE} / 1000000" | bc)

# Deleterious SNPs
# Prepare output file column names
printf '#chr\tstart_pos\tend_pos\tdel_vs_tol\tnum_del_snp\n' > ${OUT_DIR}/${OUT_PREFIX}_num_del_snps_per_${win_size_mbp}Mbp_window.txt
# Prepare windows, pull feature type, then count
bedtools makewindows -g ${GENOME_FILE} -w ${WIN_SIZE} | bedtools intersect -a - -b ${BED} -wo | awk -v Del_vs_Tol="Deleterious" '{if ($7==Del_vs_Tol) {print $0}}' | bedtools groupby -g 1,2,3,7 -c 8 -o count >> ${OUT_DIR}/${OUT_PREFIX}_num_del_snps_per_${win_size_mbp}Mbp_window.txt

# Tolerated SNPs
# Prepare output file column names
printf '#chr\tstart_pos\tend_pos\tdel_vs_tol\tnum_tol_snp\n' > ${OUT_DIR}/${OUT_PREFIX}_num_tol_snps_per_${win_size_mbp}Mbp_window.txt
# Prepare windows, pull feature type, then count
bedtools makewindows -g ${GENOME_FILE} -w ${WIN_SIZE} | bedtools intersect -a - -b ${BED} -wo | awk -v Del_vs_Tol="Tolerated" '{if ($7==Del_vs_Tol) {print $0}}' | bedtools groupby -g 1,2,3,7 -c 8 -o count >> ${OUT_DIR}/${OUT_PREFIX}_num_tol_snps_per_${win_size_mbp}Mbp_window.txt
