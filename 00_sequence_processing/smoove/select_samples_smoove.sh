#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=6gb
#SBATCH --tmp=4gb
#SBATCH -t 00:15:00
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
VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/smoove_processing/results_smoove/results_genotyped/filtered/mut_barley_cohort.filtAB-DP.noRepeatOverlap.noRefNs.vcf.gz"

# List of samples to include, must match sample names in VCF files exactly
MUT_SAMPLES="/panfs/jay/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/01_snp_filtering/Post_GATK_Filtering/wgs_mut_sample_names.txt"
HYBRID_SAMPLES="/panfs/jay/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/01_snp_filtering/Post_GATK_Filtering/wgs_hybrid_sample_names.txt"
MUT_10XGenomics_SAMPLES="/panfs/jay/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/00_sequence_processing/smoove/mut3_10xGenomics_sample_names.txt"

# Output file prefix
MUT8_PREFIX="mut8_smoove.filtAB-DP.noRepeatOverlap.noRefNs.vcf.gz"
HYBRID_PREFIX="hybrid13_smoove.filtAB-DP.noRepeatOverlap.noRefNs.vcf.gz"
MUT3_PREFIX="mut3_10xGenomics_smoove.filtAB-DP.noRepeatOverlap.noRefNs.vcf.gz"

# Output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/smoove_processing/results_smoove/results_genotyped/filtered"
SCRATCH_DIR="/scratch.global/liux1299/temp_mut8_and_hybrid_barley"

#-----------------
mkdir -p ${OUT_DIR}

# Select relevant samples and remove sites where all genotypes are missing and sites
#   that are all homozygous reference (account for missing genotypes)
bcftools view --samples-file ${MUT_SAMPLES} ${VCF} | bcftools filter -e "F_PASS(GT='mis') == 1.0 | COUNT(GT='RR')=(N_SAMPLES-COUNT(GT='mis'))" -O z -o ${OUT_DIR}/${MUT8_PREFIX}.vcf.gz
bcftools view --samples-file ${HYBRID_SAMPLES} ${VCF} | bcftools filter -e "F_PASS(GT='mis') == 1.0 | COUNT(GT='RR')=(N_SAMPLES-COUNT(GT='mis'))" -O z -o ${OUT_DIR}/${HYBRID_PREFIX}.vcf.gz
bcftools view --samples-file ${MUT_10XGenomics_SAMPLES} ${VCF} | bcftools filter -e "F_PASS(GT='mis') == 1.0 | COUNT(GT='RR')=(N_SAMPLES-COUNT(GT='mis'))" -O z -o ${OUT_DIR}/${MUT3_PREFIX}.vcf.gz

# Tabix index vcf files
parallel --verbose --tmpdir ${SCRATCH_DIR} tabix -p vcf {} ::: "${OUT_DIR}/${MUT8_PREFIX}.vcf.gz" "${OUT_DIR}/${HYBRID_PREFIX}.vcf.gz" "${OUT_DIR}/${MUT3_PREFIX}.vcf.gz"
