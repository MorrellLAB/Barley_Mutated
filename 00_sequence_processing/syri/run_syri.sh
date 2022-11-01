#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=200gb
#SBATCH --tmp=160gb
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p ram256g,ram1t,amdsmall,amdlarge,amd512,amd2tb
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

set -e
set -o pipefail

# This script utilizes Slurm job arrays, one per chromosome
#   Output files have the file prefix auto-generated as: PAF_prefix-REF.prefix-QRY.prefix.
# To run: sbatch --array=0-6 run_syri.sh
# Where "0-6" corresponds to the number of chromosomes defined in CHR_ARR (below), remember bash arrays are zero based

# Dependencies
module load python3/3.9.3_anaconda2021.11_mamba
module load gcc/7.2.0
export PATH=${PATH}:/panfs/jay/groups/9/morrellp/shared/Software/syri

# User provided input arguments
# IMPORTANT! syri depends on alignments generated with the minimap2 --eqx flag (write =/X CIGAR operators).
#   If this flag is omitted, syri will return an error similar to the following:
#   "Reading BAM/SAM file - ERROR - Incorrect CIGAR string found. CIGAR string can only have I/D/H/S/X/=."
# List of full paths to PAF files
PAF_LIST="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/morex_v2_to_v3/paf_chr_list_morex_refV3_qryV2_asm5.txt"
# List of full paths to REF FASTA files
REF_FASTA_LIST="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/assembly/split_by_chr/chr_fasta_list_HvulgareMorex_702_V3.hardmasked.txt"
# List of full paths to query FASTA files
QRY_FASTA_LIST="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/split_by_chr/chr_fasta_list_Barley_Morex_V2_pseudomolecules.txt"
# Bash array of chromosomes that match chromosome names in all three file lists above
CHR_ARR=("chr1H" "chr2H" "chr3H" "chr4H" "chr5H" "chr6H" "chr7H")
# Output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/morex_v2_to_v3/syri_out"

# Output prefix
SAMPLE="Morex"
REF_PREFIX="HvMorex_702_V3.hardmasked"
QRY_PREFIX="Morex_V2"

#----------------------
mkdir -p ${OUT_DIR}

# Slurm job array prep
# Determine maximum array limit
MAX_ARRAY_LIMIT=$[${#CHR_ARR[@]} - 1]
echo "Maximum array limit is ${MAX_ARRAY_LIMIT}."
# Get the current chromosome we are processing
CURR_CHR=${CHR_ARR[${SLURM_ARRAY_TASK_ID}]}
echo "Currently processing chromosome: ${CURR_CHR}"
# Gather set of three files for current chromosome
PAF=$(grep -w "${CURR_CHR}" ${PAF_LIST})
REF_FASTA=$(grep -w "${CURR_CHR}" ${REF_FASTA_LIST})
QRY_FASTA=$(grep -w "${CURR_CHR}" ${QRY_FASTA_LIST})
echo "Current set of files:"
echo "PAF file: ${PAF}"
echo "Reference fasta file: ${REF_FASTA}"
echo "Query fasta file: ${QRY_FASTA}"

# Generate unique output file prefix from input files
# PAF_PREFIX=$(basename ${PAF} .paf)
# # Get REF FASTA basename. Handles .fa.gz, .fasta.gz, and .fa file exensions
# if [[ ${REF_FASTA} == *".fa.gz"* ]]; then
#     REF_PREFIX=$(basename ${REF_FASTA} .fa.gz)
# elif [[ ${REF_FASTA} == *".fasta.gz"* ]]; then
#     REF_PREFIX=$(basename ${REF_FASTA} .fasta.gz)
# elif [[ "${REF_FASTA}" == *".fa" ]]; then
#     REF_PREFIX=$(basename ${REF_FASTA} .fa)
# elif [[ "${REF_FASTA}" == *".fasta" ]]; then
#     REF_PREFIX=$(basename ${REF_FASTA} .fasta)
# fi
# # Get Query FSATA basename. Handles .fa.gz, .fasta.gz, and .fa file exensions
# if [[ ${QRY_FASTA} == *".fa.gz"* ]]; then
#     QRY_PREFIX=$(basename ${QRY_FASTA} .fa.gz)
# elif [[ ${QRY_FASTA} == *".fasta.gz"* ]]; then
#     QRY_PREFIX=$(basename ${QRY_FASTA} .fasta.gz)
# elif [[ "${QRY_FASTA}" == *".fa" ]]; then
#     QRY_PREFIX=$(basename ${QRY_FASTA} .fa)
# elif [[ "${QRY_FASTA}" == *".fasta" ]]; then
#     QRY_PREFIX=$(basename ${QRY_FASTA} .fasta)
# fi
# # Auto generate output file prefix
# OUT_PREFIX=$(printf "${PAF_PREFIX}-REF.${REF_PREFIX}-QRY.${QRY_PREFIX}.")

# Build output file prefix from SAMPLE, REF_PREFIX, and QRY_PREFIX specified in user input args
OUT_PREFIX=$(printf "${CURR_CHR}_${SAMPLE}-REF.${REF_PREFIX}-QRY.${QRY_PREFIX}.")

# Run syri using PAF input format
syri -c ${PAF} -r ${REF_FASTA} -q ${QRY_FASTA} -k -F P --dir ${OUT_DIR} --prefix ${OUT_PREFIX} --lf syri_${OUT_PREFIX}log
