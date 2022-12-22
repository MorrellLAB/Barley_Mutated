#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=42gb
#SBATCH --tmp=30gb
#SBATCH -t 05:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Exclude discordant sites

# Dependencies
module load bedtools/2.29.2
module load htslib/1.9

# User provided input arguments
# Discrodant sites to exclude in sorted BED format
DISCORDANT_BED="/scratch.global/liux1299/temp_morex-sample2/temp_morex-sample2.10x_snp_indel_x_invariant.bed"
# 10x Genomics filterd phased variants SNPs and INDELs
VCF_10X_SNPs="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/diffs_from_ref/morex-sample2_phased_variants-snps.DPfilt.noRepeatOverlap.noRefNs.vcf.gz"
# Invariant sites for same sample as the 10x Genomics sample above
VCF_INVARIANT="/scratch.global/liux1299/temp_morex-sample2/morex-sample2.invariant_sites.filtDP-RGQ-QUAL.snps.vcf"
# Output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/pixy_input"
# Temporary directory
TEMP_DIR="/scratch.global/liux1299/temp_morex-sample2"

#-----------------
mkdir -p ${OUT_DIR} ${TEMP_DIR}

# Prepare output file names
BN_10X_SNPs=$(basename ${VCF_10X_SNPs} .vcf.gz)
BN_INVARIANT=$(basename ${VCF_INVARIANT} .vcf)

# Exclude discordant sites
bedtools intersect -sorted -v -header -wa -a ${VCF_10X_SNPs} -b ${DISCORDANT_BED} > ${OUT_DIR}/${BN_10X_SNPs}.concordant.vcf
bedtools intersect -sorted -v -header -wa -a ${VCF_INVARIANT} -b ${DISCORDANT_BED} > ${TEMP_DIR}/${BN_INVARIANT}.concordant.vcf

# Compress and index vcfs
bgzip ${OUT_DIR}/${BN_10X_SNPs}.concordant.vcf
tabix -p vcf ${OUT_DIR}/${BN_10X_SNPs}.concordant.vcf.gz

bgzip -c ${TEMP_DIR}/${BN_INVARIANT}.concordant.vcf > ${OUT_DIR}/${BN_INVARIANT}.concordant.vcf.gz
tabix -p vcf ${OUT_DIR}/${BN_INVARIANT}.concordant.vcf.gz
