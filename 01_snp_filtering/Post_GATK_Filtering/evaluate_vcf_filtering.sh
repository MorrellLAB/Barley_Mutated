#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=8gb
#SBATCH --tmp=6gb
#SBATCH -t 08:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

# This script uses various metrics to evaluate VCF filtering

# Dependencies
module load python3/3.8.3_anaconda2020.07_mamba
module load bcftools/1.10.2
module load texlive/20131202
export PATH=${PATH}:/panfs/jay/groups/9/morrellp/liux1299/sequence_handling/HelperScripts

# User provided input arguments
# Prior to filtering on various annotations
VCF1="/scratch.global/liux1299/temp_mut8_and_hybrid_barley/mut8_and_hybrid_barley_snps_polymorphic.vcf.gz"
# After filtering on various annotations
VCF2="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/mut8_and_hybrid_barley_snps_biallelic.vcf.gz"
VCF3="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/mut8_and_hybrid_barley_snps_biallelic.noRepeatOverlap.noRefNs.vcf.gz"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/vcf_summary"
# Arguments for bcftools stats
REF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta"

#-----------------
# Check the out dir and scratch dir exist, if not make them
mkdir -p ${OUT_DIR}

function calc_maf() {
    local vcf="$1"
    local out_dir="$2"
    # Get out_prefix
    if [[ "${vcf}" == *"gz"* ]]; then
        # We are working with gzipped VCF file
        out_prefix=$(basename ${vcf} .vcf.gz)
    else
        # We are working with uncompressed vcf
        out_prefix=$(basename ${vcf} .vcf)
    fi
    # Custom script
    echo "Starting custom script for MAF calculations for file: ${vcf}"
    VCF_MAF.py ${vcf} > ${out_dir}/${out_prefix}.maf
    echo "Done running custom MAF calculations"
}

export -f calc_maf

# bcftools stats
function generate_stats() {
    local vcf=$1
    local ref_gen=$2
    local out_dir=$3
    # Prep output file prefix
    if [[ "${vcf}" == *".gz"* ]]; then
        prefix=$(basename ${vcf} .vcf.gz)
    else
        prefix=$(basename ${vcf} .vcf)
    fi
    mkdir -p ${out_dir}/plots_${prefix}
    # Generate stats
    bcftools stats -F ${ref_gen} -s - ${vcf} > ${out_dir}/${prefix}.stats
    # Generate plots
    plot-vcfstats -p ${out_dir}/plots_${prefix} -s ${out_dir}/${prefix}.stats
}

export -f generate_stats

# Generate bcftools stats
generate_stats ${VCF3} ${REF} ${OUT_DIR}
generate_stats ${VCF2} ${REF} ${OUT_DIR}
generate_stats ${VCF1} ${REF} ${OUT_DIR}

# Calculate MAF to check filtering
calc_maf ${VCF3} ${OUT_DIR}
calc_maf ${VCF2} ${OUT_DIR}
calc_maf ${VCF1} ${OUT_DIR}
