#!/bin/bash

set -e
set -o pipefail

module load python3/3.8.3_anaconda2020.07_mamba
module load bcftools/1.10.2
module load texlive/20131202

# User provided input arguments
# Mutated lines
VCF1="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.SNPs.private.vcf.gz"
VCF2="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.INDELs.private.vcf.gz"
# hybrid lines
VCF3="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/hybrid13.SNPs.common.vcf.gz"
VCF4="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/hybrid13.INDELs.common.vcf.gz"
VCF5="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/hybrid13.SNPs.rare.vcf.gz"
VCF6="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/hybrid13.INDELs.rare.vcf.gz"
#SAMPLE_LIST="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut_sample_names.txt"
REF_GEN="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/mut_and_hybrid_bcftools_summaries"

#---------------------
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

# Check if out directory exists, if not make it
mkdir -p ${OUT_DIR}

# Generate stats for multiple cuts of the data
generate_stats ${VCF1} ${REF_GEN} ${OUT_DIR}
generate_stats ${VCF2} ${REF_GEN} ${OUT_DIR}
generate_stats ${VCF3} ${REF_GEN} ${OUT_DIR}
generate_stats ${VCF4} ${REF_GEN} ${OUT_DIR}
generate_stats ${VCF5} ${REF_GEN} ${OUT_DIR}
generate_stats ${VCF6} ${REF_GEN} ${OUT_DIR}
