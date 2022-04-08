#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=120gb
#SBATCH --tmp=100gb
#SBATCH -t 90:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p ram256g,ram1t,amd2tb
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Dependencies
module load longranger/2.2.2

# User provided input arguments
#   Unique run ID used to name project directory
SAMPLE_ID="morex-sample2"
#   Sample prefix must match FASTQ files
#   For multiple prefixes, separate with comma
#   i.e. WBDC355-1,WBDC355-2,WBDC355-3,WBDC355-4
SAMPLE_PREFIX="morex-sample2"
#   Path to the FASTQ directory
FASTQ_DIR=/panfs/roc/groups/9/morrellp/shared/Datasets/10x_Genomics/Barley/Morex
#   Full path to 10x Genomics compatible reference directory
#   Note: this is the output directory from longranger mkref
REF_DIR=
#   Specify "freebayes", path to GATK with format gatk:/paht/to/GenomeAnalysisTK.jar, or "disable"
#   Must use GATK 3.2, 3.4, or 3.5 (possibly 3.8 is ok, but we need to try to know)
#   GATK 3.6 has crashing bug in 'single-haplotype' mode
VC_MODE="gatk:/panfs/roc/groups/9/morrellp/public/Software/GATK_ML_3.8.1-0-gf15c1c3ef/GenomeAnalysisTK.jar"
#   For inbred barley, use female
#   Longranger requires this argument if it is not able to automatically detect the sex of the organism
SEX=female
#   Where do we want our output files to go?
OUT_DIR=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/${SAMPLE_ID}

#------------------------------
# Make out dir
mkdir -p ${OUT_DIR}

# Go into out dir
cd ${OUT_DIR}
# Run longranger wgs
longranger wgs --id=${SAMPLE_ID} \
               --sample=${SAMPLE_PREFIX} \
               --fastqs=${FASTQ_DIR} \
               --reference=${REF_DIR} \
               --vcmode=${VC_MODE} \
               --sex=female
