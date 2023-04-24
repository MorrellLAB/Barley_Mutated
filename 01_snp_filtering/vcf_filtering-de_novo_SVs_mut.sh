#!/bin/bash

set -e
set -o pipefail

# For larger SVs, filter down to de novo variants in mutated lines
# This script serves as a log of commands run.

# Dependencies
module load bedtools/2.29.2
module load bcftools/1.10.2
module load htslib/1.9
module load vcflib_ML/1.0.0_rc2

# User provided input arguments
VCF_DELs="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/mut_3_lines_dels_merged.callable.vcf.gz"
VCF_LARGE_SVs="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/mut_3_lines_large_svs_merged.callable.vcf.gz"
# Callable ONT INDELs
VCF_ONT_Sniffles2="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/filtered/mut_ont.private.geSup3.callable.noRefDiffs.final.DEL.vcf.gz"
VCF_ONT_cutesv="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/cutesv_calls/filtered/mut_ont_cutesv.private.callable.noRefDiffs.final.DEL.vcf.gz"

OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/de_novo_larger_svs"
#OUT_PREFIX="mut_3_lines"

# Differences from Morex reference
REF_DIFFS_10x_del="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/diffs_from_ref/morex-sample2_dels_diffs_from_ref.noBND.bed"
REF_DIFFS_ONT_DEL="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3/filtered/morex_ont.noHomRef.geSup5.callable.DEL.bed"
REF_DIFFS_ONT_INS="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3/filtered/morex_ont.noHomRef.geSup5.callable.INS.bed"
REF_DIFFS_85xONT_DEL="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_85x_ont/filtered/parts_pos/morex_85x_ont.noHomRef.geSup5.callable.parts.DEL.bed"
REF_DIFFS_85xONT_INS="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_85x_ont/filtered/parts_pos/morex_85x_ont.noHomRef.geSup5.callable.parts.INS.bed"
REF_DIFFS_PacBio_DEL="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/pacbio_morex/pacbio_morex_v3/filtered/morex_pacbio.noHomRef.geSup5.callable.DEL.bed"
REF_DIFFS_PacBio_INS="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/pacbio_morex/pacbio_morex_v3/filtered/morex_pacbio.noHomRef.geSup5.callable.INS.bed"

#---------------------
mkdir -p ${OUT_DIR} ${OUT_DIR}/Intermediates

dels_bn=$(basename ${VCF_DELs} .vcf.gz)
large_svs_bn=$(basename ${VCF_LARGE_SVs} .vcf.gz)

# Remove SVs that overlap with diffs in ref of our Morex sample and Mascher et al. Morex
# 10x Genomics dels
bedtools intersect -wa -v -header -a ${VCF_DELs} -b ${REF_DIFFS_10x_del} ${REF_DIFFS_ONT_DEL} ${REF_DIFFS_ONT_INS} ${REF_DIFFS_85xONT_DEL} ${REF_DIFFS_85xONT_INS} ${REF_DIFFS_PacBio_DEL} ${REF_DIFFS_PacBio_INS} > ${OUT_DIR}/${dels_bn}.noRefDiffs.vcf
# 10x Genomics large svs
# One of the larger svs remain after removing overlaps with diffs from ref
bedtools intersect -wa -v -header -a ${VCF_LARGE_SVs} -b ${REF_DIFFS_10x_del} ${REF_DIFFS_ONT_DEL} ${REF_DIFFS_ONT_INS} ${REF_DIFFS_85xONT_DEL} ${REF_DIFFS_85xONT_INS} ${REF_DIFFS_PacBio_DEL} ${REF_DIFFS_PacBio_INS} > ${OUT_DIR}/${large_svs_bn}.noRefDiffs.vcf

# Only include SVs private to each line
# Create a file where only one sample has an alternate allele
bcftools view -i "COUNT(GT='alt')=1" ${OUT_DIR}/${dels_bn}.noRefDiffs.vcf -O z -o ${OUT_DIR}/${dels_bn}.noRefDiffs.private.vcf.gz
tabix -p vcf ${OUT_DIR}/${dels_bn}.noRefDiffs.private.vcf.gz

# Manually curate small remaining set by scoring in SV-Plaudit
svplaudit_scored="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/samplot-mut_10x_dels/mut_10x_dels.private.callable.noMorexDiffs.supports.vcf"

# Rescue SV that got removed during SV-Plaudit scoring but has additional support from ONT
bedtools intersect -wa -a ${OUT_DIR}/${dels_bn}.noRefDiffs.private.vcf.gz -b ${VCF_ONT_Sniffles2} ${VCF_ONT_cutesv} | uniq > ${OUT_DIR}/tmp_${dels_bn}_overlap_ONT.noRefDiffs.private.vcf

# Concatenate SV-Plaudit scored 10x dels and "rescued" SV
grep "#" ${svplaudit_scored} > ${OUT_DIR}/${dels_bn}.noRefDiffs.private.supports.vcf
grep -v "#" ${svplaudit_scored} | cat - ${OUT_DIR}/tmp_${dels_bn}_overlap_ONT.noRefDiffs.private.vcf | sort -k1,1 -k2,2n | uniq >> ${OUT_DIR}/${dels_bn}.noRefDiffs.private.supports.vcf

# Add tag "BasesToClosestVariant" to remove (mostly) consecutive variants since they are likely not de novo
vcfdistance < ${OUT_DIR}/${dels_bn}.noRefDiffs.private.supports.vcf | bcftools view -e 'BasesToClosestVariant <= 10 & SVLEN < -10' -O v -o ${OUT_DIR}/${dels_bn}.noRefDiffs.private.supports.final.vcf

# Cleanup
mv ${OUT_DIR}/tmp_${dels_bn}_overlap_ONT.noRefDiffs.private.vcf \
    ${OUT_DIR}/${dels_bn}.noRefDiffs.private.supports.vcf \
    ${OUT_DIR}/*dels_merged.callable.noRefDiffs.vcf* \
    ${OUT_DIR}/Intermediates
