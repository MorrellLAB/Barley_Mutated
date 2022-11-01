#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=8gb
#SBATCH --tmp=6gb
#SBATCH -t 00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Dependencies
module load parallel/20210822
module load samtools/1.9

# User provided input arguments
# Define chr array
CHR_ARR=("chr1H" "chr2H" "chr3H" "chr4H" "chr5H" "chr6H" "chr7H")

FA_MOREX_V3="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/assembly/HvulgareMorex_702_V3.hardmasked.fa"
OUT_DIR_MOREX_V3="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/assembly/split_by_chr"

FA_MOREX_V2="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules.fasta"
OUT_DIR_MOREX_V2="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/split_by_chr"

SCRATCH_DIR="/scratch.global/liux1299/temp"

#-------------
# Make output dir if it doesn't exist
mkdir -p ${OUT_DIR_MOREX_V3} ${OUT_DIR_MOREX_V2} ${SCRATCH_DIR}

function subset_fa_by_chr() {
    local fa_file="$1"
    local curr_chr="$2"
    local out_dir="$3"
    # Get basename of fasta file
    if [[ "${fa_file}" == *".fa.gz"* ]]; then
        out_prefix=$(basename ${fa_file} .fa.gz)
    elif [[ "${fa_file}" == *".fasta.gz" ]]; then
        out_prefix=$(basename ${fa_file} .fasta.gz)
    elif [[ "${fa_file}" == *".fa" ]]; then
        out_prefix=$(basename ${fa_file} .fa)
    elif [[ "${fa_file}" == *".fasta" ]]; then
        out_prefix=$(basename ${fa_file} .fasta)
    fi
    # Subset fasta file by chromosome
    samtools faidx ${fa_file} ${curr_chr} > ${out_dir}/${out_prefix}-${curr_chr}.fa
    # Gzip fasta file
    gzip ${out_dir}/${out_prefix}-${curr_chr}.fa
}

export -f subset_fa_by_chr

# Split fasta file by chromosome
parallel --verbose --tmpdir ${SCRATCH_DIR} subset_fa_by_chr ${FA_MOREX_V3} {} ${OUT_DIR_MOREX_V3} ::: ${CHR_ARR[@]}
parallel --verbose --tmpdir ${SCRATCH_DIR} subset_fa_by_chr ${FA_MOREX_V2} {} ${OUT_DIR_MOREX_V2} ::: ${CHR_ARR[@]}
