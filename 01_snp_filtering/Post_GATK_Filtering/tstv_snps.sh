#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=22gb
#SBATCH --tmp=12gb
#SBATCH -t 05:00:00
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
# Prior to filtering on various annotations
VCF1="/scratch.global/liux1299/temp_mut8_and_hybrid_barley/mut8_and_hybrid_barley_snps_polymorphic.vcf.gz"
# After filtering on various annotations
VCF2="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/mut8_and_hybrid_barley_snps_biallelic.vcf.gz"
VCF3="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/mut8_and_hybrid_barley_snps_biallelic.noRepeatOverlap.noRefNs.vcf.gz"
# VCF shared basename for preparing output filenames
#VCF_BN="mut8_and_hybrid_barley_"
# Output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/vcf_tstv"

# Samples list, one sample per line
#MUT_SAMP="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/all_mut_samples_list.txt"
#HYBRID_SAMP="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/all_hybrid_samples_list.txt"

# Output file prefixes
# MUT_PREFIX="mut8"
# HYBRID_PREFIX="hybrid_barley"

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
# For separated mutated and hybrid lines
# vcf1_bn=$(basename ${VCF1} .vcf.gz | sed "s/${VCF_BN}//")
# vcf2_bn=$(basename ${VCF2} .vcf.gz | sed "s/${VCF_BN}//")
# vcf3_bn=$(basename ${VCF3} .vcf.gz | sed "s/${VCF_BN}//")

# Separate mutated lines from hybrid barley lines
# We'll calculate Ts/Tv for them separately
# VCF1
# bcftools view --samples-file ${MUT_SAMP} ${VCF1} -O z -o ${OUT_DIR}/${MUT_PREFIX}_${vcf1_bn}.vcf.gz
# bcftools view --samples-file ${HYBRID_SAMP} ${VCF1} -O z -o ${OUT_DIR}/${HYBRID_PREFIX}_${vcf1_bn}.vcf.gz
# # VCF2
# bcftools view --samples-file ${MUT_SAMP} ${VCF2} -O z -o ${OUT_DIR}/${MUT_PREFIX}_${vcf2_bn}.vcf.gz
# bcftools view --samples-file ${HYBRID_SAMP} ${VCF2} -O z -o ${OUT_DIR}/${HYBRID_PREFIX}_${vcf2_bn}.vcf.gz
# # VCF3
# bcftools view --samples-file ${MUT_SAMP} ${VCF3} -O z -o ${OUT_DIR}/${MUT_PREFIX}_${vcf3_bn}.vcf.gz
# bcftools view --samples-file ${HYBRID_SAMP} ${VCF3} -O z -o ${OUT_DIR}/${HYBRID_PREFIX}_${vcf3_bn}.vcf.gz

# Calculate per sample Ts/Tv
# Per sample Ts/Tv
# Redirect only summary (stdout) to output file
# VCF3 per sample
java -jar ${SNPSIFT_JAR} tstv ${VCF3} 1> ${OUT_DIR}/${vcf3_prefix}.snpsiftTsTv.txt
#java -jar ${SNPSIFT_JAR} tstv ${OUT_DIR}/${MUT_PREFIX}_${vcf3_bn}.vcf.gz 1> ${OUT_DIR}/${MUT_PREFIX}_${vcf3_bn}.snpsiftTsTv.txt
# VCF2 per sample
java -jar ${SNPSIFT_JAR} tstv ${VCF2} 1> ${OUT_DIR}/${vcf2_prefix}.snpsiftTsTv.txt
# VCF1 per sample
java -jar ${SNPSIFT_JAR} tstv ${VCF1} 1> ${OUT_DIR}/${vcf1_prefix}.snpsiftTsTv.txt

# Calculate Ts/Tv
# VCF3
# summarize_tstv ${OUT_DIR}/${MUT_PREFIX}_${vcf3_bn}.vcf.gz ${OUT_DIR}
# summarize_tstv ${OUT_DIR}/${HYBRID_PREFIX}_${vcf3_bn}.vcf.gz ${OUT_DIR}
# # VCF2
# summarize_tstv ${OUT_DIR}/${MUT_PREFIX}_${vcf2_bn}.vcf.gz ${OUT_DIR}
# summarize_tstv ${OUT_DIR}/${HYBRID_PREFIX}_${vcf2_bn}.vcf.gz ${OUT_DIR}
# # VCF1
# summarize_tstv ${OUT_DIR}/${MUT_PREFIX}_${vcf1_bn}.vcf.gz ${OUT_DIR}
# summarize_tstv ${OUT_DIR}/${HYBRID_PREFIX}_${vcf1_bn}.vcf.gz ${OUT_DIR}

# Before separating mutated and hybrid lines
parallel --verbose summarize_tstv {} ${OUT_DIR} ::: ${VCF1} ${VCF2} ${VCF3}
# summarize_tstv ${VCF3} ${OUT_DIR}
# summarize_tstv ${VCF2} ${OUT_DIR}
# summarize_tstv ${VCF1} ${OUT_DIR}