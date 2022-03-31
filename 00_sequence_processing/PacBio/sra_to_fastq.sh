#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=22gb
#SBATCH --tmp=12gb
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

# This script runs uses NCBI's SRA Toolkit to convert SRA files to FASTQ files.

# Dependencies
module load parallel/20210822
module load sratoolkit_ML/2.9.6

# User provided input arguments
SRA_LIST="/scratch.global/liux1299/sra_pacbio/sra_pacbio_list.txt"
OUT_DIR="/scratch.global/liux1299/sra_pacbio/fastq"

#--------------------
mkdir -p "${OUT_DIR}"

# Run fastq-dump
# Layout: single
parallel fastq-dump {} --outdir "${OUT_DIR}" --gzip :::: ${SRA_LIST}
