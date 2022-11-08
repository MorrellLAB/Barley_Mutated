#!/bin/bash

set -e
set -o pipefail

# Dependencies
# Load perl module
module load parallel/20210822
module load perl/modules.centos7.5.26.1
module load htslib/1.9
export PATH=$PATH:/panfs/jay/groups/9/morrellp/shared/Software/ensembl-vep-release-108.1

# User provided input arguments
# Note: VeP only works on bgzipped and tabix indexed VCF files
VCF_LIST="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/split_by_sample_SNPs_private/mut_snps_private_list.txt"
# Full filepath to GFF file
GFF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/gene_annotation/Hv_Morex.pgsb.Jul2020.sorted.parts.gff3.gz"
HC_GFF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/gene_annotation/Hv_Morex.pgsb.Jul2020.HC.sorted.parts.gff3.gz"
LC_GFF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/gene_annotation/Hv_Morex.pgsb.Jul2020.LC.sorted.parts.gff3.gz"
# Reference FASTA file
REF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta"
# What species are we running?
SPECIES="hordeum_vulgare"
# Where do we want our output files to go?
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/HC_LC_gff_SNPs_private_per_samp"
HC_OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/HC_gff_SNPs_private_per_samp"
LC_OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/LC_gff_SNPs_private_per_samp"
# Run_Vep.sh script
RUN_VEP_SCRIPT=~/GitHub/Barley_Mutated/02_analysis/VEP/Run_Vep.sh

#--------------------------------
# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR} ${HC_OUT_DIR} ${LC_OUT_DIR}

# Pull run_vep function from script
source ${RUN_VEP_SCRIPT}
# Run function
cd ${OUT_DIR}
parallel --verbose run_vep {} ${GFF} ${REF} ${SPECIES} :::: ${VCF_LIST}

cd ${HC_OUT_DIR}
parallel --verbose run_vep {} ${HC_GFF} ${REF} ${SPECIES} :::: ${VCF_LIST}

cd ${LC_OUT_DIR}
parallel --verbose run_vep {} ${LC_GFF} ${REF} ${SPECIES} :::: ${VCF_LIST}
