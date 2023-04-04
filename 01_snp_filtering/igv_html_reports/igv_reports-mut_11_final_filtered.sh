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
VCF2="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.INDELs.private.vcf.gz"
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
BAM9="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M01-3-3/outs/M01-3-3_phased_possorted_bam.bam"
BAM10="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M20-2-2/outs/M20-2-2_phased_possorted_bam.bam"
BAM11="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M29-2-2/outs/M29-2-2_phased_possorted_bam.bam"
# Only generate reports for groups of 3 sample BAM files at a time
# Define sample name groups
# SAMP_GROUP1="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/sample_names_M02_M11_M14.txt"
# SAMP_GROUP2="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/sample_names_M28_M35_M36.txt"
# SAMP_GROUP3="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/sample_names_M39_M41.txt"
# SAMP_GROUP4="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/sample_names_M01_M20_M29.txt"
# --flanking INT. Genomic region to include either side of variant; default=1000.
FLANK_SIZE="1000"
# Output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/igv_reports"
OUT_PREFIX_SNPs="mut8_and_3mut10xGenomics.SNPs.private"
OUT_PREFIX_INDELs="mut8_and_3mut10xGenomics.INDELs.private"
TEMP="/scratch.global/liux1299"

#------------------
# Create output dir
mkdir -p ${OUT_DIR} ${TEMP}

# Pull out chr1H_part1 and chr1H_part2 only so output file is smaller
bcftools view --regions "chr1H_part1" $VCF1 -O v -o ${TEMP}/chr1H_part1-${OUT_PREFIX_SNPs}.vcf
bcftools view --regions "chr1H_part1" $VCF2 -O v -o ${TEMP}/chr1H_part1-${OUT_PREFIX_INDELs}.vcf

# SNPs
create_report ${TEMP}/chr1H_part1-${OUT_PREFIX_SNPs}.vcf ${REF} --flanking ${FLANK_SIZE} --tracks ${TEMP}/chr1H_part1-${OUT_PREFIX_SNPs}.vcf ${BAM1} ${BAM2} ${BAM3} ${BAM4} ${BAM5} ${BAM6} ${BAM7} ${BAM8} ${BAM9} ${BAM10} ${BAM11} --output ${OUT_DIR}/igv_report-${OUT_PREFIX_SNPs}.html

# indels
create_report ${TEMP}/chr1H_part1-${OUT_PREFIX_INDELs}.vcf ${REF} --flanking ${FLANK_SIZE} --tracks ${TEMP}/chr1H_part1-${OUT_PREFIX_INDELs}.vcf ${BAM1} ${BAM2} ${BAM3} ${BAM4} ${BAM5} ${BAM6} ${BAM7} ${BAM8} ${BAM9} ${BAM10} ${BAM11} --output ${OUT_DIR}/igv_report-${OUT_PREFIX_INDELs}.html
