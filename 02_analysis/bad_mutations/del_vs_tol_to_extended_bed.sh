#!/bin/bash

set -e
set -o pipefail

# Convert deleterious_vs_tolerated.with_pseudo_pos.txt to extended BED format for dSNPs per codon plotting purposes. Outputs extended BED file with column 4 containing "Deleterious" vs "Tolerated" annotation.
# Prints to stdout

# Usage: ./del_vs_tol_to_extended_bed.sh [SNP_TABLE] > deleterious_vs_tolerated.bed

# Where: [SNP_TABLE] is tab delimited and contains barley pseudomolecules positions and "Deleterious" vs "Tolerated" annotations

# User provided input arguments
SNP_TABLE="$1"

# Remove header line to conform with extended BED file format
grep -v "Chr_pseudo" ${SNP_TABLE} | awk '{ FS="\t"; OFS="\t"; $32 = $31 - 1; print $30,$32,$31,$29 }' | sort -k1,1 -k2,2n
