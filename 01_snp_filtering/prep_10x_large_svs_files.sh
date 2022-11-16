#!/bin/bash

set -e
set -o pipefail

# Prepare large_svs.vcf.gz VCF so they only include SVs in large_sv_calls.bedpe
# We want to EXCLUDE low confidence candidates
# This script is meant as a log of commands run and will need to be modified if applied to other datasets

# File format descriptions:
#   large_svs.vcf.gz: Large-scale (≥30Kbp or inter-chromosomal) structural variant and CNV calls, including low confidence candidates, in VCF format
#   large_sv_calls.bedpe: Large-scale (≥30Kbp or inter-chromosomal) structural variant and CNV calls, excluding low confidence candidates, in BEDPE format
#   large_sv_candidates.bedpe: Large-scale (≥30Kbp or inter-chromosomal) structural variant and CNV calls, including low confidence candidates, in BEDPE format

# Dependencies
module load bcftools/1.10.2
module load bedtools/2.29.2

# User provided input arguments
# Morex
VCF_MOREX="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/outs/morex-sample2_large_svs.vcf.gz"
# This should be the large_sv_calls.bedpe file
BEDPE_MOREX="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/outs/morex-sample2_large_sv_calls.bedpe"
# IDs (one per line) that cause bedtools error: "Error: Invalid record in file"
# We'll exclude these. This list will need to be generated as each error pops up each attempt at the run
# Important: If no list, put "NA"
PROBLEM_CALLS_MOREX="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/outs/morex-sample2_large_svs.problem_ids_list.txt"

# M01
VCF_M01="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M01-3-3/outs/M01-3-3_large_svs.vcf.gz"
# This should be the large_sv_calls.bedpe file
BEDPE_M01="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M01-3-3/outs/M01-3-3_large_sv_calls.bedpe"
# IDs (one per line) that cause bedtools error: "Error: Invalid record in file"
# We'll exclude these. This list will need to be generated as each error pops up each attempt at the run]
# Important: If no list, put "NA"
PROBLEM_CALLS_M01="NA"

# M20
VCF_M20="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M20-2-2/outs/M20-2-2_large_svs.vcf.gz"
# This should be the large_sv_calls.bedpe file
BEDPE_M20="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M20-2-2/outs/M20-2-2_large_sv_calls.bedpe"
# IDs (one per line) that cause bedtools error: "Error: Invalid record in file"
# We'll exclude these. This list will need to be generated as each error pops up each attempt at the run]
# Important: If no list, put "NA"
PROBLEM_CALLS_M20="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M20-2-2/outs/M20-2-2_large_svs.problem_ids_list.txt"

# M29
VCF_M29="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M29-2-2/outs/M29-2-2_large_svs.vcf.gz"
# This should be the large_sv_calls.bedpe file
BEDPE_M29="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M29-2-2/outs/M29-2-2_large_sv_calls.bedpe"
# IDs (one per line) that cause bedtools error: "Error: Invalid record in file"
# We'll exclude these. This list will need to be generated as each error pops up each attempt at the run]
# Important: If no list, put "NA"
PROBLEM_CALLS_M29="NA"

# Scratch directory for intermediate files
TEMP_DIR="/scratch.global/liux1299/temp_10xGenomics"

#--------------------------
# Define functions
function prep_large_svs() {
    local vcf="$1"
    local bedpe_file="$2"
    local temp_dir="$3"
    local problem_calls="$4"
    # Prepare output file prefix and out_dir
    if [[ "${vcf}" == *"gz"* ]]; then
        # We are working with gzipped VCF file
        out_prefix_vcf=$(basename ${vcf} .vcf.gz)
    else
        # We are working with uncompressed vcf
        out_prefix_vcf=$(basename ${vcf} .vcf)
    fi
    out_prefix_bedpe=$(basename ${bedpe_file} .bedpe)
    out_dir=$(dirname ${bedpe_file})
    # Remove chrUn from bedpe file, we can't say much about these
    # Also remove DUP since it's hard to distinguish between a real dup and alignment issues
    grep -v 'chrUn\|TYPE=DUP' ${bedpe_file} > ${out_dir}/${out_prefix_bedpe}.nochrUn.noDUP.bedpe
    # Remove DUP from VCF (some positions cause "Error: Invalid record in file" from bedtools)
    # Also remove chrUn and SVTYPE=UNK (unknown SV type)
    if [[ "${problem_calls}" == "NA" ]]; then
        # No problem calls defined yet, don't include that part
        bcftools view --targets ^chrUn -e 'INFO/SVTYPE="DUP" | INFO/SVTYPE2="DUP" | INFO/SVTYPE="UNK" | INFO/SVTYPE2="UNK"' ${vcf} -O v -o ${temp_dir}/${out_prefix_vcf}.nochrUn.noDUP-UNK.vcf
    elif [[ "${problem_calls}" != "NA" ]]; then
        # Exclude problem call ids
        bcftools view --targets ^chrUn -e 'INFO/SVTYPE="DUP" | INFO/SVTYPE2="DUP" | INFO/SVTYPE="UNK" | INFO/SVTYPE2="UNK"' ${vcf} | bcftools view -e ID==@${problem_calls} -O v -o ${temp_dir}/${out_prefix_vcf}.nochrUn.noDUP-UNK.vcf
    fi
    # Prepare VCF header
    grep "#" ${temp_dir}/${out_prefix_vcf}.nochrUn.noDUP-UNK.vcf > ${out_dir}/${out_prefix_vcf}.calls.nochrUn.noDUP-UNK.vcf
    # EXCLUDE low confidence SV calls from VCF
    bedtools intersect -wa -a ${temp_dir}/${out_prefix_vcf}.nochrUn.noDUP-UNK.vcf -b ${out_dir}/${out_prefix_bedpe}.nochrUn.noDUP.bedpe | sort -uV -k1,1 -k2,2 >> ${out_dir}/${out_prefix_vcf}.calls.nochrUn.noDUP-UNK.vcf
}

export -f prep_large_svs

# Make scratch directory
mkdir -p ${TEMP_DIR}

# Prepare large SVs VCF files
prep_large_svs ${VCF_MOREX} ${BEDPE_MOREX} ${TEMP_DIR} ${PROBLEM_CALLS_MOREX}
prep_large_svs ${VCF_M01} ${BEDPE_M01} ${TEMP_DIR} ${PROBLEM_CALLS_M01}
prep_large_svs ${VCF_M20} ${BEDPE_M20} ${TEMP_DIR} ${PROBLEM_CALLS_M20}
prep_large_svs ${VCF_M29} ${BEDPE_M29} ${TEMP_DIR} ${PROBLEM_CALLS_M29}
