#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=6gb
#SBATCH --tmp=4gb
#SBATCH -t 48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

set -e
set -o pipefail

# Dependencies
module load bcftools/1.10.2
module load parallel/20210822

# User provided input arguments
vcf="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/diffs_from_ref/morex-sample2_phased_variants-snps.DPfilt.callable.vcf.gz"
# Sliding windows file generated from prep_sliding_win.sh script
# Formatted for bcftools
sliding_win_file="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/check_heterogeneity/sliding_win_3000_overlap_1000_reformatted.txt"
# For file naming, should match file in variable sliding_win_file
win_size_bp="3000"
step_size_bp="1000"
# List of chr part names
chr_list="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/chr_parts_list.intervals"
out_dir="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/check_heterogeneity"
scratch_dir="/scratch.global/liux1299/check_heterogeneity"

#-----------------
function count_hets_in_win() {
    local interval="$1"
    local vcf="$2"
    local out_dir="$3"
    local curr_chr="$4"
    local win_size_bp="$5"
    local step_size_bp="$6"
    # Count number of het in current sliding window (bp)
    nhet_in_win=$(bcftools view --regions ${interval} --no-header --genotype het ${vcf} | wc -l)
    chr=$(echo ${interval} | cut -d':' -f 1)
    win_start=$(echo ${interval} | tr ':' '\t' | tr '-' '\t' | cut -f 2)
    win_end=$(echo ${interval} | tr ':' '\t' | tr '-' '\t' | cut -f 3)
    midpoint=$(echo "(${win_end}+${win_start})/2" | bc)
    # Save to file for plotting
    #printf "${chr}\t${win_start}\t${win_end}\t${midpoint}\t${nhet_in_win}\n" > ${out_dir}/${chr}_${win_start}-${win_end}_num_het.txt
    printf "${chr}\t${win_start}\t${win_end}\t${midpoint}\t${nhet_in_win}\n" >> ${out_dir}/${curr_chr}-num_hets_sliding_win_${win_size_bp}_overlap_${step_size_bp}.txt
}

export -f count_hets_in_win

# Store chr parts names as array
chr_arr=($(cat ${chr_list}))
# Run as Slurm task array
# Get the current chr part
curr_chr="${chr_arr[${SLURM_ARRAY_TASK_ID}]}"
echo "Current chromosome we are processing is: ${curr_chr}"

# Prepare output directories and subdirectories
mkdir -p ${out_dir} #${scratch_dir} ${scratch_dir}/${curr_chr}

sw_prefix=$(basename ${sliding_win_file} .txt)
grep "${curr_chr}" ${sliding_win_file} > ${out_dir}/${curr_chr}_${sw_prefix}.txt

# Check if output file exists, if it does, remove it since we are appending
if [ -f ${out_dir}/${curr_chr}-num_hets_sliding_win_${win_size_bp}_overlap_${step_size_bp}.txt ]; then
    # File exists, remove before proceeding
    rm ${out_dir}/${curr_chr}-num_hets_sliding_win_${win_size_bp}_overlap_${step_size_bp}.txt
fi
# Calculate number of het in each window
# Can't use parallel processing that writes to individual files because downstream bash commands to
#   concatenate files will give an error, argument list too long
#parallel --verbose count_hets_in_win {} ${vcf} ${scratch_dir}/${curr_chr} :::: ${out_dir}/${curr_chr}_${sw_prefix}.txt
for i in $(cat ${out_dir}/${curr_chr}_${sw_prefix}.txt)
do
    count_hets_in_win ${i} ${vcf} ${out_dir} ${curr_chr} ${win_size_bp} ${step_size_bp}
done

# Combine separate num het files into single sorted file
# Need to use a for loop since argument list is too long for bash
#cat ${scratch_dir}/${curr_chr}/*_num_het.txt | sort -k1,1 -k2,2n > ${out_dir}/${curr_chr}-num_hets_sliding_win_${win_size_bp}_overlap_${step_size_bp}.txt
