#!/bin/bash

set -e
set -o pipefail

# Dependencies
# Load perl module
module load perl/modules.centos7.5.26.1
module load htslib/1.9
export PATH=$PATH:/panfs/jay/groups/9/morrellp/shared/Software/ensembl-vep-release-108.1

# User provided input arguments
# Note: VeP only works on bgzipped and tabix indexed VCF files
VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.SNPs.private.vcf.gz"
# Full filepath to GFF file
GFF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/gene_annotation/Hv_Morex.pgsb.Jul2020.sorted.parts.gff3.gz"
# Reference FASTA file
REF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta"
# What species are we running?
SPECIES="hordeum_vulgare"
# Where do we want our output files to go?
#OUT_DIR=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/HC_parts_gff
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/HC_LC_gff_SNPs_private_all_samples"
# Run_Vep.sh script
RUN_VEP_SCRIPT=~/GitHub/Barley_Mutated/02_analysis/VEP/Run_Vep.sh

#--------------------------------
# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR}
# Go into output directory
cd ${OUT_DIR}

# Pull run_vep function from script
source ${RUN_VEP_SCRIPT}
# Run function
run_vep ${VCF} ${GFF} ${REF} ${SPECIES}
