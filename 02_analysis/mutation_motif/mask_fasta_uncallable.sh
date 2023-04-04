#!/bin/bash

set -e
set -o pipefail

module load bedtools/2.29.2

# User provided input arguments
FASTA="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta"
UNCALLABLE_BED="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/uncallable_regions/morex_v3_combined_uncallable.nochrUn.bed"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/emboss_compseq"

#-------------------
mkdir -p ${OUT_DIR}

fasta_bn=$(basename ${FASTA} .fasta)

# Mask uncallable regions in fasta
bedtools maskfasta -fi ${FASTA} -bed ${UNCALLABLE_BED} -fo ${OUT_DIR}/${fasta_bn}.masked_uncallable.fasta
