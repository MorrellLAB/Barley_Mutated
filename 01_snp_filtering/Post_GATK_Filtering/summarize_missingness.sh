#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=12gb
#SBATCH --tmp=8gb
#SBATCH -t 03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Dependencies
module load vcftools_ML/0.1.16
module load R/4.0.4

# User provided input arguments
VCF1="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Variant_Recalibrator/mut8_and_hybrid_barley_snps.recalibrated.pass_sites.vcf.gz"
VCF2="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/mut8_and_hybrid_barley_snps_biallelic.vcf.gz"
REF="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta"
OUT_PREFIX="mut8_and_hybrid_barley_snps"
VCF1_PREFIX="RecalPass"
VCF2_PREFIX="Biallelic"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/ann_visualization"

# We are calling on some helper scripts from sequence_handling
#   and will need to define a few variables
seqhand="/panfs/roc/groups/9/morrellp/liux1299/sequence_handling"

#-----------------
# Check the out dir and scratch dir exist, if not make them
mkdir -p ${OUT_DIR} ${OUT_DIR}/Intermediates

VCF1_PREFIX=$(basename ${VCF1} .vcf.gz)

if [[ ${VCF1} == *".gz"* ]]; then
    # Use gzip flag
    vcftools --gzvcf ${VCF1} \
        --missing-indv \
        --out ${OUT_DIR}/${OUT_PREFIX}_missingness
    # Missing data per site
    vcftools --gzvcf ${VCF1} \
        --missing-site \
        --out ${OUT_DIR}/${OUT_PREFIX}_missingness
else
    # Uncompressed VCF
    vcftools --vcf ${VCF1} \
        --missing-indv \
        --out ${OUT_DIR}/${OUT_PREFIX}_missingness
    # Missing data per site
    vcftools --vcf ${VCF1} \
        --missing-site \
        --out ${OUT_DIR}/${OUT_PREFIX}_missingness
fi

# Visualize missingness
"${seqhand}/HelperScripts/graph_missingness.R" \
    ${OUT_DIR}/${VCF1_PREFIX}_missingness.imiss \
    ${OUT_DIR}/${VCF1_PREFIX}_missingness.lmiss \
    ${OUT_DIR}

if [[ ${VCF2} == *".gz"* ]]; then
    # Use gzip flag
    vcftools --gzvcf ${VCF2} \
        --missing-indv \
        --out ${OUT_DIR}/${OUT_PREFIX}_missingness
    # Missing data per site
    vcftools --gzvcf ${VCF2} \
        --missing-site \
        --out ${OUT_DIR}/${OUT_PREFIX}_missingness
else
    # Uncompressed VCF
    vcftools --vcf ${VCF2} \
        --missing-indv \
        --out ${OUT_DIR}/${OUT_PREFIX}_missingness
    # Missing data per site
    vcftools --vcf ${VCF2} \
        --missing-site \
        --out ${OUT_DIR}/${OUT_PREFIX}_missingness
fi

# Visualize missingness
"${seqhand}/HelperScripts/graph_missingness.R" \
    ${OUT_DIR}/${OUT_PREFIX}_missingness.imiss \
    ${OUT_DIR}/${OUT_PREFIX}_missingness.lmiss \
    ${OUT_DIR}
