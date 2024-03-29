#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=12gb
#SBATCH --tmp=8gb
#SBATCH -t 01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Dependencies
module load gatk/4.1.2

RAW_VCF="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Create_HC_Subset/mut8_and_hybrid_barley_raw_variants_indels.vcf"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Create_HC_Subset"
TEMP_DIR="/scratch.global/liux1299/temp"

#-----------------
mkdir -p ${OUT_DIR} ${OUT_DIR}/Intermediates ${TEMP_DIR}

vcf_filename=$(basename ${RAW_VCF} .vcf)

# Generate a table of annotations for exploration
gatk VariantsToTable \
    -V ${RAW_VCF} \
    -F CHROM -F POS -F TYPE -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -GF GQ \
    -O "${OUT_DIR}/Intermediates/${vcf_filename}.table"
