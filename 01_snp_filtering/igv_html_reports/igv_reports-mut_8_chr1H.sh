#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=56gb
#SBATCH --tmp=22gb
#SBATCH -t 03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Dependencies
module load bcftools/1.10.2
module load python3/3.8.3_anaconda2020.07_mamba
source ~/.bashrc
conda activate igvreports

# User provided input arguments
VCF1="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.SNPs.private.vcf.gz"
VCF2="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/temp_mut8_biallelic.callable.INDELs.noMorexDiffs.private.vcf.gz"
# Reference fasta
REF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta"
# BAM file(s)
BAM1="/scratch.global/liux1299/temp_mut8_and_hybrid_barley/M02-1_15-15.bam"
BAM2="/scratch.global/liux1299/temp_mut8_and_hybrid_barley/M11-3-3_22-12.bam"
BAM3="/scratch.global/liux1299/temp_mut8_and_hybrid_barley/M14-2-2_18-05.bam"
BAM4="/scratch.global/liux1299/temp_mut8_and_hybrid_barley/M28-3-1_19-21.bam"
BAM5="/scratch.global/liux1299/temp_mut8_and_hybrid_barley/M35-3-2_18-30.bam"
BAM6="/scratch.global/liux1299/temp_mut8_and_hybrid_barley/M36-1-2_03-19.bam"
BAM7="/scratch.global/liux1299/temp_mut8_and_hybrid_barley/M39-2-3_18-28.bam"
BAM8="/scratch.global/liux1299/temp_mut8_and_hybrid_barley/M41-2-1_03-23.bam"
# Output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered"
TEMP="/scratch.global/liux1299"

#------------------
# Create output dir
mkdir -p ${OUT_DIR}

vcf1_bn=$(basename $VCF1 .vcf.gz)
vcf2_bn=$(basename $VCF2 .vcf.gz)

# Pull out chr1H_part1 and chr1H_part2 only so output file is smaller
bcftools view --regions "chr1H_part1,chr1H_part2" $VCF1 -O v -o ${TEMP}/chr1H-${vcf1_bn}.vcf
bcftools view --regions "chr1H_part1,chr1H_part2" $VCF2 -O v -o ${TEMP}/chr1H-${vcf2_bn}.vcf

create_report ${TEMP}/chr1H-${vcf1_bn}.vcf ${REF} --tracks ${TEMP}/chr1H-${vcf1_bn}.vcf ${BAM1} ${BAM2} ${BAM3} --output ${OUT_DIR}/igv_report-chr1H-${vcf1_bn}.html

# indels
# Generate reports for groups of sample bam files so not too many bam files are loaded at a time
# BAM1-3
create_report ${TEMP}/chr1H-${vcf2_bn}.vcf ${REF} --tracks ${TEMP}/chr1H-${vcf2_bn}.vcf ${BAM1} ${BAM2} ${BAM3} --output ${OUT_DIR}/igv_report-chr1H-${vcf2_bn}_M02-M11-M14.html
#BAM4-6
create_report ${TEMP}/chr1H-${vcf2_bn}.vcf ${REF} --tracks ${TEMP}/chr1H-${vcf2_bn}.vcf ${BAM4} ${BAM5} ${BAM6} --output ${OUT_DIR}/igv_report-chr1H-${vcf2_bn}_M28-M35-M36.html
#BAM7-8
create_report ${TEMP}/chr1H-${vcf2_bn}.vcf ${REF} --tracks ${TEMP}/chr1H-${vcf2_bn}.vcf ${BAM7} ${BAM8} --output ${OUT_DIR}/igv_report-chr1H-${vcf2_bn}_M39-M41.html
