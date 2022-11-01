#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=20gb
#SBATCH --tmp=16gb
#SBATCH -t 00:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Dependencies
module load parallel/20210822
module load python3/3.9.3_anaconda2021.11_mamba
source activate /panfs/jay/groups/9/morrellp/liux1299/.conda/envs/plotsr_dep_env

# User provided input arguments
# List of .syri.out files
SYRI_OUT_LIST="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/morex_v2_to_v3/syri_out/syri_out_list.txt"
# REF/QRY fasta genomes text file list, one file per chromosome
GENOMES_TXT_LIST="/panfs/jay/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/00_sequence_processing/syri/genomes_plotsr_file_list.txt"
# Chromosome stored in bash array, must match chromosome naming used in SYRI_OUT_LIST and FA_GENOMES files
CHR_ARR=("chr1H" "chr2H" "chr3H" "chr4H" "chr5H" "chr6H" "chr7H")
# Output file basename
OUT_BASENAME="Morex_v2_vs_v3"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/morex_v2_to_v3/Morex_plotsr"
TMP_DIR="/scratch.global/liux1299/morex_v2_to_v3/Morex_plotsr/temp"

#---------------
mkdir -p ${OUT_DIR} ${TMP_DIR}

# Build lines for --sr flag (one for each .syri.out file)
#sr_lines=$(sed -e 's,^,--sr ,' ${SYRI_OUT_LIST} | tr '\n' ' ')

function run_plotsr() {
    local curr_chr="$1"
    local syri_out_list="$2"
    local genomes_txt_list="$3"
    local out_basename="$4"
    local out_dir="$5"
    # Get files for current chromosome
    syri_out=$(grep "${curr_chr}" ${syri_out_list})
    fa_genomes=$(grep "${curr_chr}" ${genomes_txt_list})
    echo "Current chromosome: ${curr_chr}"
    echo "Current syri.out file: ${syri_out}"
    echo "Current fasta genomes txt file: ${fa_genomes}"
    # Run plotsr
    plotsr \
        --sr ${syri_out} \
        --genomes ${fa_genomes} \
        -o ${out_dir}/${curr_chr}.${out_basename}.png
}

export -f run_plotsr

parallel --verbose --tmpdir ${TMP_DIR} run_plotsr {} ${SYRI_OUT_LIST} ${GENOMES_TXT_LIST} ${OUT_BASENAME} ${OUT_DIR} ::: ${CHR_ARR[@]}
