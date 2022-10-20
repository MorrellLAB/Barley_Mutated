#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=6gb
#SBATCH --tmp=4gb
#SBATCH -t 00:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Dependencies
module load bcftools/1.10.2
module load htslib/1.9
module load parallel/20210822

# Separate mutated and hybrid parents samples and select hybrid parents relevant to this study
VCF_SNP="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/mut8_and_hybrid_barley_snps_biallelic.noRepeatOverlap.noRefNs.vcf.gz"
VCF_INDEL="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/mut8_and_hybrid_barley_indels_noComplex.noRepeatOverlap.noRefNs.vcf.gz"

# List of samples to include, must match sample names in VCF files exactly
MUT_SAMPLES="/panfs/jay/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/01_snp_filtering/Post_GATK_Filtering/wgs_mut_sample_names.txt"
HYBRID_SAMPLES="/panfs/jay/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/01_snp_filtering/Post_GATK_Filtering/wgs_hybrid_sample_names.txt"

# Output file prefix
MUT_SNP_PREFIX="mut8_snps_biallelic.noRepeatOverlap.noRefNs"
HYBRID_SNP_PREFIX="hybrid13_snps_biallelic.noRepeatOverlap.noRefNs"
MUT_INDEL_PREFIX="mut8_indels_noComplex.noRepeatOverlap.noRefNs"
HYBRID_INDEL_PREFIX="hybrid13_indels_noComplex.noRepeatOverlap.noRefNs"

# Output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered"
SCRATCH_DIR="/scratch.global/liux1299/temp_mut8_and_hybrid_barley"

#----------------
mkdir -p ${OUT_DIR}

# Select relevant samples and remove sites where all genotypes are missing and sites
#   that are all homozygous reference (account for missing genotypes)
# SNPs
bcftools view --samples-file ${MUT_SAMPLES} ${VCF_SNP} | bcftools filter -e "F_PASS(GT='mis') == 1.0 | COUNT(GT='RR')=(N_SAMPLES-COUNT(GT='mis'))" -O z -o ${OUT_DIR}/${MUT_SNP_PREFIX}.vcf.gz
bcftools view --samples-file ${HYBRID_SAMPLES} ${VCF_SNP} | bcftools filter -e "F_PASS(GT='mis') == 1.0 | COUNT(GT='RR')=(N_SAMPLES-COUNT(GT='mis'))" -O z -o ${OUT_DIR}/${HYBRID_SNP_PREFIX}.vcf.gz

# Indels
bcftools view --samples-file ${MUT_SAMPLES} ${VCF_INDEL} | bcftools filter -e "F_PASS(GT='mis') == 1.0 | COUNT(GT='RR')=(N_SAMPLES-COUNT(GT='mis'))" -O z -o ${OUT_DIR}/${MUT_INDEL_PREFIX}.vcf.gz
bcftools view --samples-file ${HYBRID_SAMPLES} ${VCF_INDEL} | bcftools filter -e "F_PASS(GT='mis') == 1.0 | COUNT(GT='RR')=(N_SAMPLES-COUNT(GT='mis'))" -O z -o ${OUT_DIR}/${HYBRID_INDEL_PREFIX}.vcf.gz

# Tabix index vcf files
parallel --verbose --tmpdir ${SCRATCH_DIR} tabix -p vcf {} ::: "${OUT_DIR}/${MUT_SNP_PREFIX}.vcf.gz" "${OUT_DIR}/${HYBRID_SNP_PREFIX}.vcf.gz" "${OUT_DIR}/${MUT_INDEL_PREFIX}.vcf.gz" "${OUT_DIR}/${HYBRID_INDEL_PREFIX}.vcf.gz"
