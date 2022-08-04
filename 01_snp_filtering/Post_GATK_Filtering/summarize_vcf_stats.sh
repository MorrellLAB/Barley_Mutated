#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=8gb
#SBATCH --tmp=6gb
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
module load python3/3.6.3_anaconda5.0.1 # This version works with bcftools, v3.7.3 results in errors
module load R/3.6.3

# User provided input arguments
VCF="/scratch.global/liux1299/temp_mut8_and_hybrid_barley/mut8_and_hybrid_barley_snps_polymorphic.vcf.gz"
REF="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta"
OUT_PREFIX="mut8_and_hybrid_barley_snps"
# Scratch directory used to store intermediate files that don't need to be kept long term
SCRATCH_DIR="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/vcf_summary"
seqhand="/panfs/roc/groups/9/morrellp/liux1299/sequence_handling"

#-----------------
# Check the out dir and scratch dir exist, if not make them
mkdir -p ${SCRATCH_DIR}

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
    # This causes an error but the plots still get output
    #   Error: Neither pdflatex or tectonic were found in your PATH...
    # So, we'll leave this commented out for now.
    #plot-vcfstats -p ${out_dir}/plots_${prefix} -s ${out_dir}/${prefix}.stats
}

export -f generate_stats

# Generate VCF stats
generate_stats ${VCF} ${REF} ${SCRATCH_DIR}

#   Create a file with minor allele frequencies (MAF)
python3 "${seqhand}/HelperScripts/VCF_MAF.py" "${VCF}" > "${SCRATCH_DIR}/${OUT_PREFIX}_MAF.txt"
#   Use R to plot the MAF file as a histogram
Rscript "${seqhand}/HelperScripts/plot_maf.R" "${SCRATCH_DIR}/${OUT_PREFIX}_MAF.txt" "${OUT_PREFIX}" "${SCRATCH_DIR}/${OUT_PREFIX}_MAF.pdf"
