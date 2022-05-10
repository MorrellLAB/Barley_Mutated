#!/bin/bash

set -e
set -o pipefail

# Prepare VCF files for Samplot. Reheader and merge VCF files.
# This script serves as a log of commands run

# Dependencies
module load bcftools/1.10.2
module load htslib/1.9
module load bedtools/2.29.2

# User provided input arguments
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration/samplot_vcfs"
OUT_NAME="morex_combined"
REF="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules_parts.fasta"
REF_Ns_BED="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration/Barley_Morex_V2_pseudomolecules_parts_missing.bed"
# VCF files, one per variable in the same order as the BAM files
VCF1="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/morex-sample2/Filtered/morex-sample2_filtered_pass1.vcf"
VCF2="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_nanopore_V2_partsRef/Morex_nanopore_V2_partsRef_90.vcf"
VCF3="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/pacbio_morex/pacbio_morex_v2/Morex_pacbio_90.vcf"
# New sample names
NEW_NAME_10x="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration/samplot_vcfs/new_name_10x.txt"
#NEW_NAME_ONT="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration/samplot_vcfs/new_name_ont.txt"
NEW_NAME_PacBio="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/exploration/samplot_vcfs/new_name_pacbio.txt"

#----------------
# Make output directories
mkdir -p ${OUT_DIR}

# Rename sample names for each VCF
# 10x Genomics VCF
OUT_PREFIX1=$(basename ${VCF1} .vcf)
bcftools reheader -s ${NEW_NAME_10x} ${VCF1} > ${OUT_DIR}/${OUT_PREFIX1}_renamed.vcf
# Nanopore VCF
#OUT_PREFIX2=$(basename ${VCF2} .vcf)
#bcftools reheader -s ${NEW_NAME_ONT} ${VCF2} > ${OUT_DIR}/${OUT_PREFIX2}_renamed.vcf
# PacBio VCF
OUT_PREFIX3=$(basename ${VCF3} .vcf)
bcftools reheader -s ${NEW_NAME_PacBio} ${VCF3} > ${OUT_DIR}/${OUT_PREFIX3}_renamed.vcf

# Exclude regions with long stretches of Ns in the REF
# 10x Genomics VCF
grep "#" ${OUT_DIR}/${OUT_PREFIX1}_renamed.vcf > ${OUT_DIR}/${OUT_PREFIX1}_renamed_noMissing.vcf
bedtools intersect -a ${OUT_DIR}/${OUT_PREFIX1}_renamed.vcf -b ${REF_Ns_BED} -v >> ${OUT_DIR}/${OUT_PREFIX1}_renamed_noMissing.vcf
# Nanopore VCF
grep "#" ${OUT_DIR}/${OUT_PREFIX2}_renamed.vcf > ${OUT_DIR}/${OUT_PREFIX2}_renamed_noMissing.vcf
bedtools intersect -a ${OUT_DIR}/${OUT_PREFIX2}_renamed.vcf -b ${REF_Ns_BED} -v >> ${OUT_DIR}/${OUT_PREFIX2}_renamed_noMissing.vcf
# PacBio VCF
grep "#" ${OUT_DIR}/${OUT_PREFIX3}_renamed.vcf > ${OUT_DIR}/${OUT_PREFIX3}_renamed_noMissing.vcf
bedtools intersect -a ${OUT_DIR}/${OUT_PREFIX3}_renamed.vcf -b ${REF_Ns_BED} -v >> ${OUT_DIR}/${OUT_PREFIX3}_renamed_noMissing.vcf

# Gzip and index VCFs
bgzip ${OUT_DIR}/${OUT_PREFIX1}_renamed_noMissing.vcf
bgzip ${OUT_DIR}/${OUT_PREFIX2}_renamed_noMissing.vcf
bgzip ${OUT_DIR}/${OUT_PREFIX3}_renamed_noMissing.vcf

VCF_10x="${OUT_DIR}/${OUT_PREFIX1}_renamed_noMissing.vcf.gz"
VCF_ONT="${OUT_DIR}/${OUT_PREFIX2}_renamed_noMissing.vcf.gz"
VCF_PacBio="${OUT_DIR}/${OUT_PREFIX3}_renamed_noMissing.vcf.gz"

tabix -p vcf ${VCF_10x}
tabix -p vcf ${VCF_ONT}
tabix -p vcf ${VCF_PacBio}

# Make sure there aren't issues with REF alleles mismatching
# 10x Genomics VCF
bcftools norm -c w -f ${REF} -o ${OUT_DIR}/temp_ref_check_${OUT_PREFIX1}.vcf -O v ${VCF_10x} >& temp_ref_check_warn_${OUT_PREFIX1}.log
# Nanopore VCF
bcftools norm -c w -f ${REF} -o ${OUT_DIR}/temp_ref_check_${OUT_PREFIX2}.vcf -O v ${VCF_ONT} >& temp_ref_check_warn_${OUT_PREFIX2}.log
# PacBio VCF
bcftools norm -c w -f ${REF} -o ${OUT_DIR}/temp_ref_check_${OUT_PREFIX3}.vcf -O v ${VCF_PacBio} >& temp_ref_check_warn_${OUT_PREFIX3}.log

# Merge VCF files
bcftools merge ${VCF_10x} ${VCF_ONT} ${VCF_PacBio} -o ${OUT_DIR}/${OUT_NAME}.vcf -O v --threads 6
