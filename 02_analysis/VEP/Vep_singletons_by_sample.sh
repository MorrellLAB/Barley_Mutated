#!/bin/bash

set -e
set -o pipefail

# Dependencies
# Load perl module
module load perl/modules.centos7.5.26.1
module load parallel
export PATH=$PATH:/panfs/roc/groups/9/morrellp/liux1299/Software/ensembl-vep

# User provided input arguments
# Note: VeP only works on bgzipped and tabix indexed VCF files
VCF_LIST=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated/Filtered/by_sample/mut_3_lines_gz_vcf_list.txt
# Full filepath to GFF file
GFF=/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_gene_annotation_PGSB.all.parts.gff3.gz
# Reference FASTA file
REF=/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules_parts.fasta
# What species are we running?
SPECIES=hordeum_vulgare
# Where do we want our output files to go?
<<<<<<< HEAD
OUT_DIR=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/per_sample-singletons
=======
OUT_DIR=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/per_sample
>>>>>>> 693ad1ebff45606877ec5674c076a4ad7a6c87bf

#--------------------------------
# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR}
# Go into output directory
cd ${OUT_DIR}

# Prepare vcf file array
VCF_ARR=($(cat ${VCF_LIST}))

function run_vep() {
    local vcf=$1
    local gff_file=$2
    local ref=$3
    local species=$4
    # Prepare out prefix from filename
    out_prefix=$(basename ${vcf} .vcf.gz)
    # Run VeP
    vep -i ${vcf} \
        --format vcf \
        --check_svs \
        --custom ${gff_file},,gff \
        --fasta ${ref} \
        --gff ${gff_file} \
        --output_file ${out_prefix} \
        --species ${species} \
        --total_length \
        --verbose \
        --warning_file ${out_prefix}.log
    # Rename txt file output
    mv ${out_prefix} ${out_prefix}.txt
}

export -f run_vep

parallel run_vep {} ${GFF} ${REF} ${SPECIES} ::: ${VCF_ARR[@]}
