#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=56gb
#SBATCH --tmp=40gb
#SBATCH -t 42:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

# Use HiFiAdapterFilt tool to filter adapters from PacBio data downloaded
#   from NCBI SRA:
#   Software: https://github.com/sheinasim/HiFiAdapterFilt
#   Paper: Sim et al. 2022 BMC Genomics. https://doi.org/10.1186/s12864-022-08375-1

set -e
set -o pipefail

# Dependencies
module load bamtools/2.5.1
export PATH=${PATH}:/panfs/roc/groups/9/morrellp/liux1299/Software/ncbi-blast-2.13.0+/bin
export PATH=${PATH}:/panfs/roc/groups/9/morrellp/liux1299/Software/HiFiAdapterFilt-2.0.0
export PATH=${PATH}:/panfs/roc/groups/9/morrellp/liux1299/Software/HiFiAdapterFilt-2.0.0/DB

# User provided input arguments
# Directory containing FASTQ files
FASTQ_DIR="/scratch.global/liux1299/sra_pacbio/fastq"
OUT_DIR="/scratch.global/liux1299/sra_pacbio/adapter_filtered_fastq"
NUM_THREADS="16"

#--------------------
# Make output dir
mkdir -p ${OUT_DIR}
# Go into directory containing all fastq files to run
cd ${FASTQ_DIR}

# Run HiFiAdapterFilt on all fastq files in ${FASTQ_DIR}
bash pbadapterfilt.sh -t ${NUM_THREADS} -o ${OUT_DIR}
