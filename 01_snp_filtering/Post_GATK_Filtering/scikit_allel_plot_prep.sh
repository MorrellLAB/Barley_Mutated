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
unfiltered_vcf="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/scikit_allel_files/chr1H_mut8_and_hybrid_barley_snps_polymorphic.vcf"
known_snp_vcf="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/scikit_allel_files/known_x_chr1H_mut8_and_hybrid_barley_snps_polymorphic.vcf"
filt1_snp_vcf="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/scikit_allel_files/chr1H_mut8_and_hybrid_barley_snps_biallelic.vcf"
filt2_snp_vcf="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/scikit_allel_files/chr1H_mut8_and_hybrid_barley_snps_biallelic.noRepeatOverlap.noRefNs.vcf"
work_dir="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/scikit_allel_files"
out_prefix_known_ann="ann_known-chr1H_mut8_and_hybrid_barley_snps_polymorphic"
out_prefix_filt_ann="ann_filt-chr1H_mut8_and_hybrid_barley_snps_polymorphic"
out_prefix_final="ann_known_and_filt-chr1H_mut8_and_hybrid_barley_snps_polymorphic"

# Randomly sample variants for plotting purposes
chrom_prefix="chr1H"
ref_fai="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta.fai"
gen_num="1000000"
gen_len="100"

#---------------------
# Go into working directory
cd ${work_dir}

# Create new tab delimited file with new annotation added
# First, add known vs novel variant annotations
vcf_to_bcftools_ann_tsv-known.py ${unfiltered_vcf} ${known_snp_vcf} > ${work_dir}/${out_prefix_known_ann}.txt

# Next, add filtered vs retained_filt1 vs retained_filt2 annotations
vcf_to_bcftools_ann_tsv-filtered.py ${unfiltered_vcf} ${filt1_snp_vcf} ${filt2_snp_vcf} > ${work_dir}/${out_prefix_filt_ann}.txt

# bgzip and tabix index before running bcftools annotate
bgzip ${work_dir}/${out_prefix_known_ann}.txt
tabix -p vcf ${work_dir}/${out_prefix_known_ann}.txt.gz
bgzip ${work_dir}/${out_prefix_filt_ann}.txt
tabix -p vcf ${work_dir}/${out_prefix_filt_ann}.txt.gz

# Prepare .hdr files containing VCF header lines to be added to output vcf
echo '##INFO=<ID=VAR_KNOWN,Number=1,Type=String,Description="Custom variant category, novel vs known">' > ${work_dir}/${out_prefix_known_ann}.hdr
echo '##INFO=<ID=VAR_FILT,Number=1,Type=String,Description="Custom variant category corresponding to each level of filtering, filtered vs retained_filt1 vs retained_filt2">' > ${work_dir}/${out_prefix_filt_ann}.hdr

# Add new TAG for known vs novel variant to INFO field
bcftools annotate -a ${work_dir}/${out_prefix_known_ann}.txt.gz -h ${work_dir}/${out_prefix_known_ann}.hdr -c CHROM,POS,REF,ALT,-,VAR_KNOWN ${unfiltered_vcf} > ${work_dir}/${out_prefix_known_ann}.vcf

# Add filtered TAG for filtered vs retained_filt1 vs retained_filt2 to INFO field
bcftools annotate -a ${work_dir}/${out_prefix_filt_ann}.txt.gz -h ${work_dir}/${out_prefix_filt_ann}.hdr -c CHROM,POS,REF,ALT,-,VAR_FILT ${work_dir}/${out_prefix_known_ann}.vcf > ${work_dir}/${out_prefix_final}.vcf

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
bedtools intersect -header -wa -a ${work_dir}/${out_prefix_final}.vcf -b ${work_dir}/${chrom_prefix}_Genome_Random_Intervals_${suffix}.bed > ${work_dir}/${chrom_prefix}_Genome_Random_Intervals_${suffix}.vcf

# Convert VCF to hdf5 format
vcf_to_h5.py ${work_dir}/${out_prefix_final}.vcf
vcf_to_h5.py ${work_dir}/${chrom_prefix}_Genome_Random_Intervals_${suffix}.vcf

# Remove intermediate file with duplicate info
rm ${work_dir}/${out_prefix_known_ann}.vcf
