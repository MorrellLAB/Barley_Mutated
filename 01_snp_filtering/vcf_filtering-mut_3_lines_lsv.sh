#!/bin/bash

set -e
set -o pipefail

# Filter 10x Genomics VCF file containing 3 mutated lines (includes SNPs, dels, and SVs)
# Output VCF file containing singletons

# Dependencies
module load bcftools/1.9
module load htslib/1.9
module load bedtools/2.29.2
module load vcftools_ML/0.1.16
module load gatk/4.1.2
# Required for vcf-annotate that is part of vcftools
export PERL5LIB=$PERL5LIB:/panfs/jay/groups/9/morrellp/public/Software/vcftools_ML-0.1.16/share/perl5

# User provided input arguments
# M01
DEL_VCF_M01="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M01-3-3/outs/M01-3-3_dels.vcf.gz"
LSV_VCF_M01="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/fix_ref_mismatch/M01-3-3_large_svs_noRefMismatch.vcf.gz"
# M20
DEL_VCF_M20="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M20-2-2/outs/M20-2-2_dels.vcf.gz"
LSV_VCF_M20="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/fix_ref_mismatch/M20-2-2_large_svs_noRefMismatch.vcf.gz"
# M29
DEL_VCF_M29="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M29-2-2/outs/M29-2-2_dels.vcf.gz"
LSV_VCF_M29="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/fix_ref_mismatch/M29-2-2_large_svs_noRefMismatch.vcf.gz"
# Full filepath to output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered"
# Output file prefix
PREFIX_M01="M01-3-3"
PREFIX_M20="M20-2-2"
PREFIX_M29="M29-2-2"
PREFIX="mut_3_lines"

# Allele balance filter, minimum and maximum cutoff
# AB filter applied to phased variants only because on this VCF has the AD format field
MIN_AB="0.30"
MAX_AB="0.70"

VCFTOOLS_CUSTOM_FILTER="~/GitHub/Barley_Mutated/01_snp_filtering/filters.txt"

# Uncallable BED file includes: REF has stretches of N's, repeat annotations, and high copy regions
UNCALLABLE_BED="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/uncallable_regions/morex_v3_combined_uncallable.nochrUn.bed"
# Uncallable BED file that includes the above plus low complexity regions
UNCALLABLE_LOW_COMPLEXITY_BED="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/uncallable_regions/morex_v3_combined_uncallable.low_complexity.nochrUn.bed"
# High diversity, >2% diversity in a 400bp window for morex-sample2
HIGH_DIV_BED="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/uncallable_regions/pixy_pi_400bp_win.gt0.02.bed"

# Differences from Morex reference
# BED file containing sites that differ between 10x morex-sample2 and Morex reference
REF_DIFFS_10x_del="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/diffs_from_ref/morex-sample2_dels_diffs_from_ref.noBND.bed"
REF_DIFFS_ONT_DEL="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3/filtered/morex_ont.noHomRef.geSup5.callable.DEL.bed"
REF_DIFFS_ONT_INS="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_ont_partsRefv3/filtered/morex_ont.noHomRef.geSup5.callable.INS.bed"
REF_DIFFS_85xONT_DEL="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_85x_ont/filtered/parts_pos/morex_85x_ont.noHomRef.geSup5.callable.parts.DEL.bed"
REF_DIFFS_85xONT_INS="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_85x_ont/filtered/parts_pos/morex_85x_ont.noHomRef.geSup5.callable.parts.INS.bed"
REF_DIFFS_PacBio_DEL="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/pacbio_morex/pacbio_morex_v3/filtered/morex_pacbio.noHomRef.geSup5.callable.DEL.bed"
REF_DIFFS_PacBio_INS="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/pacbio_morex/pacbio_morex_v3/filtered/morex_pacbio.noHomRef.geSup5.callable.INS.bed"

#----------------------------------------
# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR} ${OUT_DIR}/Intermediates-mut

# Prepare large SVs VCF before filtering
function prep_large_svs() {
    local large_svs_vcf="$1"
    # Prepare large SVs VCF basename and dirname
    lsv_dir=$(dirname ${large_svs_vcf})
    if [[ "${large_svs_vcf}" == *"vcf.gz" ]]; then
        lsv_bn=$(basename ${large_svs_vcf} .vcf.gz)
    elif [[ "${large_svs_vcf}" == *"vcf" ]]; then
        lsv_bn=$(basename ${large_svs_vcf} .vcf)
    fi
    # These will then be stored in bash arrays
    echo "${lsv_dir}"
    echo "${lsv_bn}"
    # Separate: INV, DEL, and BND
    # These will have different filtering criteria
    bcftools view -i 'INFO/SVTYPE="DEL"' ${large_svs_vcf} -O z -o ${lsv_dir}/${lsv_bn}.DEL.vcf.gz
    bcftools view -i 'INFO/SVTYPE="INV"' ${large_svs_vcf} -O z -o ${lsv_dir}/${lsv_bn}.INV.vcf.gz
    bcftools view -i 'INFO/SVTYPE="BND" | INFO/SVTYPE2="BND"' ${large_svs_vcf} -O z -o ${lsv_dir}/${lsv_bn}.BND.vcf.gz
    # Index vcfs
    tabix -p vcf ${lsv_dir}/${lsv_bn}.DEL.vcf.gz
    tabix -p vcf ${lsv_dir}/${lsv_bn}.INV.vcf.gz
    tabix -p vcf ${lsv_dir}/${lsv_bn}.BND.vcf.gz
}

export -f prep_large_svs

# Prepare large SVs before additional filtering
LSV_M01_ARR=($(prep_large_svs ${LSV_VCF_M01}))
LSV_M01_DIR=${LSV_M01_ARR[0]}
LSV_M01_BN=${LSV_M01_ARR[1]}

LSV_M20_ARR=($(prep_large_svs ${LSV_VCF_M20}))
LSV_M20_DIR=${LSV_M20_ARR[0]}
LSV_M20_BN=${LSV_M20_ARR[1]}

LSV_M29_ARR=($(prep_large_svs ${LSV_VCF_M29}))
LSV_M29_DIR=${LSV_M29_ARR[0]}
LSV_M29_BN=${LSV_M29_ARR[1]}

# Merge samples for each type of variant
# deletions VCF
bcftools merge ${DEL_VCF_M01} ${DEL_VCF_M20} ${DEL_VCF_M29} -O z -o ${OUT_DIR}/Intermediates-mut/${PREFIX}_dels.vcf.gz
# larger SVs - deletions
bcftools merge ${LSV_M01_DIR}/${LSV_M01_BN}.DEL.vcf.gz ${LSV_M20_DIR}/${LSV_M20_BN}.DEL.vcf.gz ${LSV_M29_DIR}/${LSV_M29_BN}.DEL.vcf.gz -O z -o ${OUT_DIR}/Intermediates-mut/${PREFIX}_larger_svs.DEL.vcf.gz
# larger SVs - inversions
bcftools merge ${LSV_M01_DIR}/${LSV_M01_BN}.INV.vcf.gz ${LSV_M20_DIR}/${LSV_M20_BN}.INV.vcf.gz ${LSV_M29_DIR}/${LSV_M29_BN}.INV.vcf.gz -O z -o ${OUT_DIR}/Intermediates-mut/${PREFIX}_larger_svs.INV.vcf.gz

# Filter out sites using 10x Genomics custom filters and exclude chrUn
# and separate BND variants
# deletions VCF
# no BND
bcftools filter --targets "^chrUn" -e 'FILTER=="LOWQ"' ${OUT_DIR}/Intermediates-mut/${PREFIX}_dels.vcf.gz | bcftools view -e 'INFO/SVTYPE=="BND" | INFO/SVTYPE2=="BND"' -O z -o ${OUT_DIR}/Intermediates-mut/${PREFIX}_dels.10xCustomFilt.noBND.vcf.gz
# BND only
bcftools filter --targets "^chrUn" -e 'FILTER=="LOWQ"' ${OUT_DIR}/Intermediates-mut/${PREFIX}_dels.vcf.gz | bcftools view -i 'INFO/SVTYPE=="BND" | INFO/SVTYPE2=="BND"' -O z -o ${OUT_DIR}/Intermediates-mut/${PREFIX}_dels.10xCustomFilt.BND.vcf.gz
# Large SVs vcf
# DEL
bcftools filter --targets "^chrUn" -e 'FILTER=="LOWQ"' ${OUT_DIR}/Intermediates-mut/${PREFIX}_larger_svs.DEL.vcf.gz | bcftools filter -e 'ALT~"chrUn"' -O z -o ${OUT_DIR}/Intermediates-mut/${PREFIX}_large_svs.DEL.10xCustomFilt.vcf.gz
# INV
bcftools filter --targets "^chrUn" -e 'FILTER=="LOWQ"' ${OUT_DIR}/Intermediates-mut/${PREFIX}_larger_svs.INV.vcf.gz | bcftools filter -e 'ALT~"chrUn"' -O z -o ${OUT_DIR}/Intermediates-mut/${PREFIX}_large_svs.INV.10xCustomFilt.vcf.gz

# Only keep variants private (i.e. present in only one sample) to each sample
# deletions VCF
bcftools view -i "COUNT(GT='alt')=1" ${OUT_DIR}/Intermediates-mut/${PREFIX}_dels.10xCustomFilt.noBND.vcf.gz -O z -o ${OUT_DIR}/Intermediates-mut/${PREFIX}_dels.10xCustomFilt.noBND.private.vcf.gz
# Larger SVs vcf
# DEL
bcftools view -i "COUNT(GT='alt')=1" ${OUT_DIR}/Intermediates-mut/${PREFIX}_large_svs.DEL.10xCustomFilt.vcf.gz -O z -o ${OUT_DIR}/Intermediates-mut/${PREFIX}_large_svs.DEL.10xCustomFilt.private.vcf.gz
# INV
bcftools view -i "COUNT(GT='alt')=1" ${OUT_DIR}/Intermediates-mut/${PREFIX}_large_svs.INV.10xCustomFilt.vcf.gz -O z -o ${OUT_DIR}/Intermediates-mut/${PREFIX}_large_svs.INV.10xCustomFilt.private.vcf.gz

# Remove SVs that overlap with uncallable regions and morex differences
# deletions VCF
bedtools intersect -wa -v -header -a ${OUT_DIR}/Intermediates-mut/${PREFIX}_dels.10xCustomFilt.noBND.private.vcf.gz -b ${UNCALLABLE_LOW_COMPLEXITY_BED} ${REF_DIFFS_10x_del} ${REF_DIFFS_ONT_DEL} ${REF_DIFFS_85xONT_DEL} ${REF_DIFFS_PacBio_DEL} | bgzip > ${OUT_DIR}/${PREFIX}_dels.private.callable.noMorexDiffs.vcf.gz
# Larger SVs vcf
# DEL
bedtools intersect -wa -v -header -a ${OUT_DIR}/Intermediates-mut/${PREFIX}_large_svs.DEL.10xCustomFilt.private.vcf.gz -b ${UNCALLABLE_BED} ${REF_DIFFS_10x_del} ${REF_DIFFS_ONT_DEL} ${REF_DIFFS_85xONT_DEL} ${REF_DIFFS_PacBio_DEL} | bgzip > ${OUT_DIR}/${PREFIX}_large_svs.DEL.private.callable.noMorexDiffs.vcf.gz
# INV
# Note: the 2 inversions in morex-sample2 (after filtering) overlap with uncallable regions, so no INV were identified that are diffs from ref in morex
# So, we just use uncallable here since there are 0 morex-sample2 INV diffs
bedtools intersect -wa -v -header -a ${OUT_DIR}/Intermediates-mut/${PREFIX}_large_svs.INV.10xCustomFilt.private.vcf.gz -b ${UNCALLABLE_BED} | bgzip > ${OUT_DIR}/${PREFIX}_large_svs.INV.private.callable.noMorexDiffs.vcf.gz
# Check variants with igv-reports (see igv_html_reports subdirectory for commands/scripts)
