#!/bin/bash

set -e
set -o pipefail

# Dependencies
# Load perl module
module load perl/modules.centos7.5.26.1
export PATH=$PATH:/panfs/roc/groups/9/morrellp/liux1299/Software/ensembl-vep

# User provided input arguments
# Note: VeP only works on bgzipped and tabix indexed VCF files
VCF=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated/Filtered/mut_3_lines_filtered_hom_singletons_only.vcf.gz
# Output file prefix
OUT_PREFIX=Morex_Mutants-hom_singletons_only
# Full filepath to GFF file
GFF=/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/gene_annotation/Barley_Morex_V2_gene_annotation_PGSB.HC.parts.gff3.gz
# Reference FASTA file
REF=/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules_parts.fasta
# What species are we running?
SPECIES=hordeum_vulgare
# Where do we want our output files to go?
OUT_DIR=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/HC_parts_gff

#--------------------------------
# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR}
# Go into output directory
cd ${OUT_DIR}

# Run VeP
vep -i ${VCF} \
    --format vcf \
    --check_svs \
    --custom ${GFF},,gff \
    --fasta ${REF} \
    --gff ${GFF} \
    --output_file ${OUT_PREFIX} \
    --species ${SPECIES} \
    --total_length \
    --verbose \
    --warning_file ${OUT_PREFIX}.log

# Rename txt file output
mv ${OUT_PREFIX} ${OUT_PREFIX}.txt
