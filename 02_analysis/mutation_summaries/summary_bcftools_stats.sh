#!/bin/bash

set -e
set -o pipefail

module load bcftools/1.10.2
module load python3/3.6.3_anaconda5.0.1
module load texlive/20131202

# User provided input arguments
VCF1="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated/Filtered/mut_3_lines_filtered_singletons_only.vcf.gz"
VCF2="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated/Filtered/mut_3_lines_filtered_no_morex_diffs.vcf"
VCF3="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/morex-sample2/Filtered/morex-sample2_filtered_pass2_for_bcftools.vcf"
VCF_BY_SAMPLE_LIST="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated/Filtered/by_sample/mut_3_lines_gz_vcf_list.txt"
REF_GEN="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules_parts.fasta"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/mutation_summaries/bcftools_summaries"

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

function generate_stats_by_sample() {
    local vcf=$1
    local ref_gen=$2
    local out_dir=$3
    # Prepare sample array
    curr_sample=$(zgrep "#CHROM" ${vcf} | cut -f 10)
    echo "Currently processing sample: ${curr_sample}"
    # Prep output file prefix
    if [[ "${vcf}" == *".gz"* ]]; then
        prefix=$(basename ${vcf} .vcf.gz)
    else
        prefix=$(basename ${vcf} .vcf)
    fi
    mkdir -p ${out_dir}/plots_${prefix}
    # Generate stats
    bcftools stats -F ${ref_gen} -s ${curr_sample} ${vcf} > ${out_dir}/${prefix}.stats
    # Generate plots
    plot-vcfstats -p ${out_dir}/plots_${prefix} -s -t ${curr_sample} -T ${curr_sample} ${out_dir}/${prefix}.stats
}

export -f generate_stats_by_sample

# Check if out directory exists, if not make it
mkdir -p ${OUT_DIR}

# Generate stats for multiple cuts of the data
generate_stats ${VCF1} ${REF_GEN} ${OUT_DIR}
generate_stats ${VCF2} ${REF_GEN} ${OUT_DIR}
generate_stats ${VCF3} ${REF_GEN} ${OUT_DIR}

# Generate stats by samples
# Mutated lines
VCF_ARR=($(cat ${VCF_BY_SAMPLE_LIST}))
for i in ${VCF_ARR[@]}
do
    generate_stats_by_sample ${i} ${REF_GEN} ${OUT_DIR}
done
