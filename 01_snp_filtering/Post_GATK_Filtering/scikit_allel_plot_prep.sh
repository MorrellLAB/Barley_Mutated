#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=22gb
#SBATCH --tmp=12gb
#SBATCH -t 01:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# This script prepares a tab delimited .txt file containing the new annotation to add
#   to the INFO field. The output .txt file is formatted to work with bcftools annotate
# The resulting output files make it easier for variant exploration using scikit allel and
#   python's plotting tools

# Dependencies
module load htslib/1.9
module load bcftools/1.10.2
module load bedtools/2.29.2
# Make sure cyvcf2 is installed for current python version
module load python3/3.8.3_anaconda2020.07_mamba
export PATH=${PATH}:"/panfs/jay/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/01_snp_filtering/Post_GATK_Filtering"
export PATH=${PATH}:"/panfs/jay/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/01_snp_filtering/vcf_exploration"

# User provided input arguments
# VCF with custom annotations subset to chr1H only
vcf_ann="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/scikit_allel_files/ann_known_and_filt-chr1H_mut8_and_hybrid_barley_snps_polymorphic.vcf"
vcf_ann_known="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/scikit_allel_files/only_known-mut8_and_hybrid_barley_snps_polymorphic.vcf"
work_dir="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/scikit_allel_files"

# Randomly sample variants for plotting purposes
chrom_prefix="chr1H"
ref_fai="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta.fai"
gen_num="1000000"
gen_len="100"

#---------------------
# Go into working directory
cd ${work_dir}

# Downsample variants (to speed up computational time)
# Create .bed file with random genome intervals
awk -v OFS='\t' {'print $1,$2'} "${ref_fai}" | grep "${chrom_prefix}" > ${work_dir}/${chrom_prefix}_genome_file.txt
# Assign variables for filenames
num_name=$(echo "scale=2; $gen_num/1000000" | bc)
if [[ "$(echo $num_name | cut -c 2-)" == ".00" ]]; then
    num_name=$(echo "scale=0; $gen_num/1000000" | bc)
fi
suffix="${num_name}Mx${gen_len}bp"
# Create intervals file to subset genome randomly
bedtools random -l ${gen_len} -n ${gen_num} -seed 65 -g "${work_dir}/${chrom_prefix}_genome_file.txt" | sort -k 1,1 -k2,2n > "${work_dir}/${chrom_prefix}_Genome_Random_Intervals_${suffix}.bed"
# Create VCF from random intervals file
bedtools intersect -header -wa -a ${vcf_ann} -b ${work_dir}/${chrom_prefix}_Genome_Random_Intervals_${suffix}.bed > ${work_dir}/${chrom_prefix}_Genome_Random_Intervals_${suffix}.vcf

# Convert VCF to hdf5 format
vcf_to_h5.py ${vcf_ann}
vcf_to_h5.py ${work_dir}/${chrom_prefix}_Genome_Random_Intervals_${suffix}.vcf
vcf_to_h5.py ${vcf_ann_known}
