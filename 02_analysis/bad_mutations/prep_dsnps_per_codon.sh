#!/bin/bash

set -e
set -o pipefail

# Run all steps needed to prepare files for generating dsnps per codon plot
# Note: this is specific to barley parts positions, for other datasets, this will need modifications

# Dependencies
module load python3/3.8.3_anaconda2020.07_mamba
# Path to directory containing scripts called below
export PATH=${PATH}:~/GitHub/Barley_Mutated/02_analysis/bad_mutations

# User provided input arguments
# deleterious_vs_tolerated.txt file parts positions
DEL_VS_TOL="$1"
# Fasta index (.fai) file works, pseudomolecules positions
GENOME_FILE="$2"
# Positions of transcripts from GFF pseudomolecules positions
GFF="$3"
# Window size in bp
WIN_SIZE="$4"
OUT_DIR="$5"
OUT_PREFIX="$6"

#----------
mkdir -p ${OUT_DIR}

# Convert Deleterious vs Tolerated file to pseudomolecules positions
del_vs_tol_basedir=$(dirname ${DEL_VS_TOL})
del_vs_tol_bn=$(basename ${DEL_VS_TOL} .txt)

# Convert parts to pseudomolecules positions
del_vs_tol_barley_parts_to_pseudomolecules.py ${DEL_VS_TOL} > ${del_vs_tol_basedir}/${del_vs_tol_bn}.with_pseudo_pos.txt
# Convert pseudo position and classification to BED format
del_vs_tol_to_extended_bed.sh ${del_vs_tol_basedir}/${del_vs_tol_bn}.with_pseudo_pos.txt > ${del_vs_tol_basedir}/${del_vs_tol_bn}.bed

# Count number of "Deleterious" and "Tolerated" SNPs per 10 Mb window
num_del_vs_tol_per_window.sh ${GENOME_FILE} ${del_vs_tol_basedir}/${del_vs_tol_bn}.bed ${OUT_DIR} ${OUT_PREFIX} ${WIN_SIZE}
