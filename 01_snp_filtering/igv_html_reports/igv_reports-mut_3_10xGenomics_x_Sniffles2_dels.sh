#!/bin/bash

set -e
set -o pipefail

# Dependencies
module load bcftools/1.10.2
module load python3/3.8.3_anaconda2020.07_mamba
source ~/.bashrc
conda activate igvreports

# User provided input arguments
VCF1="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/de_novo_larger_svs/ont_support/concat_mut_3_lines_dels_10xGenomics_and_ONTSniffles2.callable.noRefDiffs.final.vcf"
# Reference fasta
REF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta"
# BAM file(s)
BAM1="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M01-3-3/outs/M01-3-3_phased_possorted_bam.bam"
BAM2="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M20-2-2/outs/M20-2-2_phased_possorted_bam.bam"
BAM3="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M29-2-2/outs/M29-2-2_phased_possorted_bam.bam"
BAM4="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/M01_ont_partsRefv3/M01_ont_partsRefv3_90_wRG.bam"
BAM5="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/M20_ont_partsRefv3/M20_ont_partsRefv3_90_wRG.bam"
BAM6="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/M29_ont_partsRefv3/M29_ont_partsRefv3_90_wRG.bam"
# --flanking INT. Genomic region to include either side of variant; default=1000.
FLANK_SIZE="2000"
# Output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/igv_reports"
OUT_PREFIX="concat_mut_3_lines_dels_10xGenomics_and_ONTSniffles2.callable.noRefDiffs.final"
TEMP="/scratch.global/liux1299"

#------------------
# Create output dir
mkdir -p ${OUT_DIR}

create_report ${VCF1} ${REF} --flanking ${FLANK_SIZE} --tracks ${VCF1} ${BAM1} ${BAM2} ${BAM3} ${BAM4} ${BAM5} ${BAM6} --output ${OUT_DIR}/igv_report-${OUT_PREFIX}.html
