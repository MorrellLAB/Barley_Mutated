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
VCF1="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/mut_3_lines_phased_variants.private.callable.noMorexDiffs.SNPs.vcf.gz"
VCF2="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/mut_3_lines_phased_variants.private.callable.noMorexDiffs.indels.vcf.gz"
# Reference fasta
REF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta"
# BAM file(s)
BAM1="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M01-3-3/outs/M01-3-3_phased_possorted_bam.bam"
BAM2="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M20-2-2/outs/M20-2-2_phased_possorted_bam.bam"
BAM3="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M29-2-2/outs/M29-2-2_phased_possorted_bam.bam"
# Output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered"
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

create_report ${TEMP}/chr1H-${vcf2_bn}.vcf ${REF} --tracks ${TEMP}/chr1H-${vcf2_bn}.vcf ${BAM1} ${BAM2} ${BAM3} --output ${OUT_DIR}/igv_report-chr1H-${vcf2_bn}.html
