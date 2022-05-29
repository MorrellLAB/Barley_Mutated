#!/bin/bash

# This script serves as a log of commands run

set -e
set -o pipefail

# Dependencies
module load bcftools/1.9

# User provided input arguments
# List of VCF files to merge
# Filtered Morex 10x Genomics vcf
VCF_10x="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/morex-sample2_all_var_filt_concat.vcf.gz"
# Filtered Morex ONT vcf
VCF_ONT="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3/filtered/morex_ont_noHomRef_geSup3.vcf.gz"
# Filtered Morex pacbio vcf
VCF_PACBIO="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/pacbio_morex/pacbio_morex_v3/filtered/morex_pacbio_noHomRef_geSup3.vcf.gz"
# Output directory
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration_morex_v3"
# Output file prefix
PREFIX="morex_10x_ont_pacbio_noHomRef"

#---------------------
# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR} ${OUT_DIR}/ref_allele_differs

# Merge VCF files
# There may be cases where the REF alleles differ at the same position.
# Example: chr1H_part2:230606101 was called as a:
#   SNP in the 10x Genomics VCF
#   Insertion in the ONT VCF
#   Didn't exist in the PacBio VCF
# We'll deal with these on a case by case basis by removing them from the
#   merged VCF and saving them to a separate file to investigate separately
# Exclude sites where REF alleles differ
printf "chr1H_part2\t230606101\n" > ${OUT_DIR}/ref_allele_differs/ref_allele_differs_regions_list.txt
printf "chr6H_part2\t189228325\n" >> ${OUT_DIR}/ref_allele_differs/ref_allele_differs_regions_list.txt

# Prepare 10x VCF with excluded REF alleles differ sites
BASENAME_10x=$(basename ${VCF_10x} .vcf.gz)
bcftools view --targets-file ^"${OUT_DIR}/ref_allele_differs/ref_allele_differs_regions_list.txt" ${VCF_10x} -O z -o ${OUT_DIR}/ref_allele_differs/${BASENAME_10x}_noRefMismatch.vcf.gz
# Index vcf
bcftools index ${OUT_DIR}/ref_allele_differs/${BASENAME_10x}_noRefMismatch.vcf.gz
# Store file in variable
vcf_10x_clean="${OUT_DIR}/ref_allele_differs/${BASENAME_10x}_noRefMismatch.vcf.gz"
# Pull sites where REF alleles differ into separate files for exploration later on
bcftools view --targets-file "${OUT_DIR}/ref_allele_differs/ref_allele_differs_regions_list.txt" ${VCF_10x} -O v -o ${OUT_DIR}/ref_allele_differs/${BASENAME_10x}_RefAlleleMismatch.vcf

# Prepare ONT vcf with excluded REF alleles differ sites
BASENAME_ONT=$(basename ${VCF_ONT} .vcf.gz)
bcftools view --targets-file ^"${OUT_DIR}/ref_allele_differs/ref_allele_differs_regions_list.txt" ${VCF_ONT} -O z -o ${OUT_DIR}/ref_allele_differs/${BASENAME_ONT}_noRefMismatch.vcf.gz
# Index vcf
bcftools index ${OUT_DIR}/ref_allele_differs/${BASENAME_ONT}_noRefMismatch.vcf.gz
# Store file in variable
vcf_ont_clean="${OUT_DIR}/ref_allele_differs/${BASENAME_ONT}_noRefMismatch.vcf.gz"
# Pull sites where REF alleles differ into separate files for exploration later on
bcftools view --targets-file "${OUT_DIR}/ref_allele_differs/ref_allele_differs_regions_list.txt" ${VCF_ONT} -O v -o ${OUT_DIR}/ref_allele_differs/${BASENAME_ONT}_RefAlleleMismatch.vcf

# Merge VCF files
bcftools merge -m both ${vcf_10x_clean} ${vcf_ont_clean} ${VCF_PACBIO} --threads 6 | bcftools sort -O z -o ${OUT_DIR}/${PREFIX}.vcf.gz
# Index VCF
bcftools index ${OUT_DIR}/${PREFIX}.vcf.gz
