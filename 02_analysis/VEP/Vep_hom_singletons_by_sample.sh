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
VCF_LIST=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated/Filtered/by_sample/mut_3_lines_hom_singletons_gz_vcf_list.txt
# Full filepath to GFF file
GFF=/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/gene_annotation/Barley_Morex_V2_gene_annotation_PGSB.HC.parts.gff3.gz
# Reference FASTA file
REF=/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules_parts.fasta
# What species are we running?
SPECIES=hordeum_vulgare
# Where do we want our output files to go?
OUT_DIR=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/HC_parts_gff/per_sample-hom_singletons
# Run_Vep.sh script
RUN_VEP_SCRIPT=~/GitHub/Barley_Mutated/02_analysis/VEP/Run_Vep.sh

#--------------------------------
# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR}
# Go into output directory
cd ${OUT_DIR}

# Prepare vcf file array
VCF_ARR=($(cat ${VCF_LIST}))

# Pull run_vep function from script
source ${RUN_VEP_SCRIPT}
# Run function
parallel run_vep {} ${GFF} ${REF} ${SPECIES} ::: ${VCF_ARR[@]}
