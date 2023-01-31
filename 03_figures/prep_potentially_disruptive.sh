#!/bin/bash

set -e
set -o pipefail

# Prepare files necessary for plotting

# Dependencies
module load bedtools/2.29.2

# User provided input arguments
# GFF Morex v3 parts positions
GFF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/annotation/HvulgareMorex_702_V3.gene.parts.nochrUn.gff3"
# VCFs
MUT_SNPs="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.SNPs.private.vcf.gz"
MUT_INDELs="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.INDELs.private.vcf.gz"
HYB_RARE_SNPs="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/hybrid13.SNPs.rare.vcf.gz"
HYB_RARE_INDELs="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/hybrid13.INDELs.rare.vcf.gz"
HYB_COMM_SNPs="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/hybrid13.SNPs.common.vcf.gz"
HYB_COMM_INDELs="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/hybrid13.INDELs.common.vcf.gz"

#---------------
function prep_file() {
    local vcf=$1
    local gff=$2
    # Get working dir
    work_dir=$(dirname ${vcf})
    if [[ ${vcf} == *".gz" ]]; then
        out_prefix=$(basename ${vcf} .vcf.gz)
    else
        out_prefix=$(basename ${vcf} .vcf)
    fi
    # Get variants that overlap with GFF
    bedtools intersect -header -wa -a ${vcf} -b ${gff} | uniq > ${work_dir}/${out_prefix}.gffOverlap.vcf
}

export -f prep_file

prep_file ${MUT_SNPs} ${GFF}
prep_file ${MUT_INDELs} ${GFF}

prep_file ${HYB_RARE_SNPs} ${GFF}
prep_file ${HYB_RARE_INDELs} ${GFF}

prep_file ${HYB_COMM_SNPs} ${GFF}
prep_file ${HYB_COMM_INDELs} ${GFF}
