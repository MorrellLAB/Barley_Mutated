#!/bin/bash

set -e
set -o pipefail

# Dependencies
module load bcftools/1.10.2
module load python3/3.8.3_anaconda2020.07_mamba
source ~/.bashrc
conda activate igvreports

# User provided input arguments
VCF1="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3/filtered/morex_ont.noHomRef.geSup5.callable.INDELs.vcf"
# Reference fasta
REF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta"
# BAM file(s)
BAM1="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3/Morex_ont_partsRefv3_90_wRG.bam"
# --flanking INT. Genomic region to include either side of variant; default=1000.
FLANK_SIZE="2000"
# Output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/igv_reports"
OUT_PREFIX="sniffles2-morex_ont.noHomRef.geSup5.callable.INDELs"
TEMP="/scratch.global/liux1299"

#------------------
# Create output dir
mkdir -p ${OUT_DIR}

create_report ${VCF1} ${REF} --flanking ${FLANK_SIZE} --tracks ${VCF1} ${BAM1} --output ${OUT_DIR}/igv_report-${OUT_PREFIX}.html
