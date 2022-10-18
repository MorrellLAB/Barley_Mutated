#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=12gb
#SBATCH --tmp=8gb
#SBATCH -t 02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Dependencies
module load gatk/4.1.2

RAW_VCF="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Create_HC_Subset/mut8_and_hybrid_barley_raw_variants.vcf"
HC_SUBSET_VCF="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Create_HC_Subset/mut8_and_hybrid_barley_high_confidence_subset.vcf"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Create_HC_Subset"
TEMP_DIR="/scratch.global/liux1299/temp"

#-----------------
mkdir -p ${OUT_DIR} ${TEMP_DIR}

vcf_filename=$(basename ${RAW_VCF} .vcf)
hc_vcf_filename=$(basename ${HC_SUBSET_VCF} .vcf)

# Select indels only
gatk SelectVariants \
    -V ${RAW_VCF} \
    -select-type INDEL \
    -O "${OUT_DIR}/${vcf_filename}_indels.vcf" \
    --tmp-dir ${TEMP_DIR}

# Select indels only from high confidence subset
gatk SelectVariants \
    -V ${HC_SUBSET_VCF} \
    -select-type INDEL \
    -O "${OUT_DIR}/${hc_vcf_filename}_indels.vcf" \
    --tmp-dir ${TEMP_DIR}
