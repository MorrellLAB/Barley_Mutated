#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=22gb
#SBATCH --tmp=12gb
#SBATCH -t 01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

# This script serves as a log for calculating Ts/Tv of sodium azide
#   lines separately from hybrid barley lines.

# Dependencies
module load parallel/20210822
module load vcftools_ML/0.1.16
module load java/openjdk-17.0.2
SNPSIFT_JAR="/panfs/jay/groups/9/morrellp/shared/Software/snpEff/SnpSift.jar"

# User provided input arguments
# 3 10xGenomics mutated lines and 8 WGS mutated lines, SNPs private to each line
VCF1="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.SNPs.private.vcf.gz"
# Hybrid WGS lines
VCF2="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/hybrid13_snps_biallelic.callable.vcf.gz"
# 3 10x Genomics mutated lines and 8 WGS mutated lines
VCF3="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/Intermediates/mut8_and_3mut10xGenomics.SNPs.vcf.gz"
# Hybrid WGS lines - rare
VCF4="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/hybrid13.SNPs.rare.vcf.gz"

# Output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/tstv_final_filtered"

#-----------------
# Check the out dir and scratch dir exist, if not make them
mkdir -p ${OUT_DIR}

function summarize_tstv() {
    local vcf="$1"
    local scratch_dir="$2"
    if [[ ${vcf} == *".gz"* ]]; then
        out_prefix=$(basename ${vcf} .vcf.gz)
        # Use gzip flag
        vcftools --gzvcf ${vcf} --FILTER-summary --out ${scratch_dir}/${out_prefix}
        vcftools --gzvcf ${vcf} --TsTv-summary --out ${scratch_dir}/${out_prefix}
    else
        out_prefix=$(basename ${vcf} .vcf)
        vcftools --vcf ${vcf} --FILTER-summary --out ${scratch_dir}/${out_prefix}
        vcftools --vcf ${vcf} --TsTv-summary --out ${scratch_dir}/${out_prefix}
    fi
}

export -f summarize_tstv

# Prepare output filenames
vcf1_prefix=$(basename ${VCF1} .vcf.gz)
vcf2_prefix=$(basename ${VCF2} .vcf.gz)
vcf3_prefix=$(basename ${VCF3} .vcf.gz)
vcf4_prefix=$(basename ${VCF4} .vcf.gz)

# Calculate per sample Ts/Tv
# Per sample Ts/Tv
# Redirect only summary (stdout) to output file
# VCF1 per sample
java -jar ${SNPSIFT_JAR} tstv ${VCF1} 1> ${OUT_DIR}/${vcf1_prefix}.snpsiftTsTv.txt
# VCF2 per sample
java -jar ${SNPSIFT_JAR} tstv ${VCF2} 1> ${OUT_DIR}/${vcf2_prefix}.snpsiftTsTv.txt
# VCF3 per sample
java -jar ${SNPSIFT_JAR} tstv ${VCF3} 1> ${OUT_DIR}/${vcf3_prefix}.snpsiftTsTv.txt
# VCF4 per sample
java -jar ${SNPSIFT_JAR} tstv ${VCF4} 1> ${OUT_DIR}/${vcf4_prefix}.snpsiftTsTv.txt

# Calculate Ts/Tv overall summary
parallel --verbose summarize_tstv {} ${OUT_DIR} ::: ${VCF1} ${VCF2} ${VCF3} ${VCF4}
