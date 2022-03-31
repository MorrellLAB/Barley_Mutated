#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=180gb
#SBATCH --tmp=90gb
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

# This script runs a quality assessment on PacBio raw FASTQ files and
#   checks for adapter contamination. Submits as job arrays, one sample per array.
# Usage: sbatch --array=0-4 quality_assessment_raw_fastq.sh

set -e
set -o pipefail

# Dependencies
module load python3/3.7.4_anaconda2019.10
source activate /panfs/roc/groups/9/morrellp/liux1299/.conda/envs/longqc_env

# User provided input arguments
FASTQ_LIST="/scratch.global/liux1299/sra_pacbio/fastq/Morex_pacbio_fastq_list.txt"
PLATFORM_KIT="pb-sequel"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/pacbio_morex"
# Dir that contains longQC.py script and modified minimap2
LONGQC_DIR="/panfs/roc/groups/9/morrellp/liux1299/Software/LongQC"

#--------------------
# Make output dir
mkdir -p ${OUT_DIR}
# Go into dir containing longQC.py script
cd ${LONGQC_DIR}

# Store current list of samples in array
FASTQ_ARR=($(cat ${FASTQ_LIST}))
# Determine maximum array limit
MAX_ARRAY_LIMIT=$[${#FASTQ_ARR[@]} - 1]
echo "Maximum array limit is ${MAX_ARRAY_LIMIT}."
# Get the current sample we are processing
CURR_FASTQ=${FASTQ_ARR[${SLURM_ARRAY_TASK_ID}]}

function run_longqc_sample_qc() {
    local fastq_file="$1"
    local platform_kit="$2"
    local out_dir="$3"
    echo "Current fastq file: ${fastq_file}"
    out_prefix=$(basename ${fastq_file} .1.fastq.gz)
    # Run Sample QC
    python longQC.py sampleqc -x ${platform_kit} \
        -o ${out_dir}/${out_prefix} \
        ${fastq_file}
}

export -f run_longqc_sample_qc

run_longqc_sample_qc ${CURR_FASTQ} ${PLATFORM_KIT} ${OUT_DIR}
#parallel --verbose run_longqc_sample_qc {} ${PLATFORM_KIT} ${OUT_DIR} :::: ${FASTQ_LIST}
