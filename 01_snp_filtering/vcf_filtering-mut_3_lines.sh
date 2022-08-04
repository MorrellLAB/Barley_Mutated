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
# Required for vcf-annotate that is part of vcftools
export PERL5LIB=$PERL5LIB:/panfs/roc/groups/9/morrellp/public/Software/vcftools_ML-0.1.16/share/perl5

# User provided input arguments
#RAW_VCF=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated/mut_3_lines_sorted.vcf
# M01
DEL_VCF_M01="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M01-3-3/outs/M01-3-3_dels.vcf.gz"
LSV_VCF_M01="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/fix_ref_mismatch/M01-3-3_large_svs_noRefMismatch.vcf.gz"
PHV_VCF_M01="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/complex_variants/M01-3-3_phased_variants_noComplex.vcf.gz"
# M20
DEL_VCF_M20="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M20-2-2/outs/M20-2-2_dels.vcf.gz"
LSV_VCF_M20="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/fix_ref_mismatch/M20-2-2_large_svs_noRefMismatch.vcf.gz"
PHV_VCF_M20="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/complex_variants/M20-2-2_phased_variants_noComplex.vcf.gz"
# M29
DEL_VCF_M29="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M29-2-2/outs/M29-2-2_dels.vcf.gz"
LSV_VCF_M29="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/fix_ref_mismatch/M29-2-2_large_svs_noRefMismatch.vcf.gz"
PHV_VCF_M29="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/complex_variants/M29-2-2_phased_variants_noComplex.vcf.gz"
# Full filepath to output directory
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered"
# Output file prefix
PREFIX_M01="M01-3-3"
PREFIX_M20="M20-2-2"
PREFIX_M29="M29-2-2"
PREFIX="mut_3_lines"

VCFTOOLS_CUSTOM_FILTER="~/GitHub/Barley_Mutated/01_snp_filtering/filters.txt"

# List of regions where REF has stretches of N's
REF_Ns_BED="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/stretches_of_Ns/Barley_MorexV3_pseudomolecules_parts_missing.bed"
# Repeat annotations
REPEAT_ANN="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/annotation/HvulgareMorex_702_V3.repeatmasked_assembly_V3.parts.gff3"
# High copy regions (e.g., chloroplasts, mitochondria, rDNA repeats, centromere repeats, etc.)
HIGH_COPY_BED="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/high_copy_regions/Morex_v3_high_copy_uniq.parts.bed"

# BED file containing sites that differ between 10x Morex and Morex reference
#BED_EXCLUSION_LIST="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/morex-sample2/Filtered/morex-sample2_diffs_from_ref.bed"

#----------------------------------------
# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR}

# Pass 1 filtering
# Filter out sites using 10x Genomics custom filters and exclude chrUn
function pass1_filtering() {
    local del_vcf="$1"
    local large_svs_vcf="$2"
    local phased_var_vcf="$3"
    local out_dir="$4"
    local sample_prefix="$5"
    ### DEL
    bcftools filter --targets "^chrUn" -e 'FILTER=="LOWQ"' ${del_vcf} -O z -o ${out_dir}/${sample_prefix}_dels.10xCustomFilt.vcf.gz
    # Separate BNDs from other variant types
    # File of other SVs, excluded BNDs
    bcftools view -e 'INFO/SVTYPE=="BND" | INFO/SVTYPE2=="BND"' ${out_dir}/${sample_prefix}_dels.10xCustomFilt.vcf.gz -O z -o ${out_dir}/${sample_prefix}_dels.10xCustomFilt.noBND.vcf.gz
    # File of only BND variants
    bcftools view -i 'INFO/SVTYPE=="BND" | INFO/SVTYPE2=="BND"' ${out_dir}/${sample_prefix}_dels.10xCustomFilt.vcf.gz -O z -o ${out_dir}/${sample_prefix}_dels.10xCustomFilt.BND_only.vcf.gz
    ### Large SVs
    bcftools filter --targets "^chrUn" -e 'FILTER=="LOWQ"' ${large_svs_vcf} | bcftools filter -e 'ALT~"chrUn"' -O z -o ${out_dir}/${sample_prefix}_large_svs.10xCustomFilt.vcf.gz
    # Separate BNDs from other variant types
    # File of other SVs, excluded BNDs
    bcftools view -e 'INFO/SVTYPE=="BND" | INFO/SVTYPE2=="BND" | INFO/SVTYPE=="UNK"' ${out_dir}/${sample_prefix}_large_svs.10xCustomFilt.vcf.gz -O z -o ${out_dir}/${sample_prefix}_large_svs.10xCustomFilt.noBND.vcf.gz
    # File of only BND variants
    bcftools view -i 'INFO/SVTYPE=="BND" | INFO/SVTYPE2=="BND" | INFO/SVTYPE=="UNK"' ${out_dir}/${sample_prefix}_large_svs.10xCustomFilt.vcf.gz -O z -o ${out_dir}/${sample_prefix}_large_svs.10xCustomFilt.BND_only.vcf.gz
    # Phased variants
    # Use 10x Genomics tuned filters and filter out sites with DP < 5
    bcftools filter --targets "^chrUn" -e 'FILTER=="10X_QUAL_FILTER" || FILTER=="10X_ALLELE_FRACTION_FILTER" || FILTER=="10X_PHASING_INCONSISTENT" || FILTER=="10X_HOMOPOLYMER_UNPHASED_INSERTION" || FILTER=="10X_RESCUED_MOLECULE_HIGH_DIVERSITY"' ${phased_var_vcf} | bcftools filter -e 'INFO/DP < 5' -O z -o ${out_dir}/${sample_prefix}_phased_variants.10xCustomFilt.vcf.gz
}

export -f pass1_filtering

function pass2_filtering() {
    local del_vcf="$1"
    local large_svs_vcf="$2"
    local phased_var_vcf="$3"
    local out_dir="$4"
    local sample_prefix="$5"
    local repeat_ann="$6"
    local high_copy_bed="$7"
    local ref_ns_bed="$8"
    # Remove SVs that overlap with repeat annotated regions, high copy regions, and that overlap with stretches of Ns
    # DEL
    bedtools intersect -wa -v -header -a ${del_vcf} -b ${repeat_ann} ${high_copy_bed} ${ref_ns_bed} | bgzip > ${out_dir}/${sample_prefix}_dels.noRepeatOverlap.noRefNs.vcf.gz
    # Large SVs
    bedtools intersect -wa -v -header -a ${large_svs_vcf} -b ${repeat_ann} ${high_copy_bed} ${ref_ns_bed} | bgzip > ${out_dir}/${sample_prefix}_large_svs.noRepeatOverlap.noRefNs.vcf.gz
    # Phased variants
    bedtools intersect -wa -v -header -a ${phased_var_vcf} -b ${repeat_ann} ${high_copy_bed} ${ref_ns_bed} | bgzip > ${out_dir}/${sample_prefix}_phased_variants.noRepeatOverlap.noRefNs.vcf.gz
    # Index vcfs
    tabix -p vcf ${out_dir}/${sample_prefix}_dels.noRepeatOverlap.noRefNs.vcf.gz
    tabix -p vcf ${out_dir}/${sample_prefix}_large_svs.noRepeatOverlap.noRefNs.vcf.gz
    tabix -p vcf ${out_dir}/${sample_prefix}_phased_variants.noRepeatOverlap.noRefNs.vcf.gz
}

export -f pass2_filtering

# Pass1 filtering: 10x Genomics custom filters and separating BND variants
pass1_filtering ${DEL_VCF_M01} ${LSV_VCF_M01} ${PHV_VCF_M01} ${OUT_DIR} ${PREFIX_M01}
pass1_filtering ${DEL_VCF_M20} ${LSV_VCF_M20} ${PHV_VCF_M20} ${OUT_DIR} ${PREFIX_M20}
pass1_filtering ${DEL_VCF_M29} ${LSV_VCF_M29} ${PHV_VCF_M29} ${OUT_DIR} ${PREFIX_M29}

# Pass2 filtering
# Remove SVs that overlap with repeat annotated regions and that overlap with stretches of Ns
# M01
pass2_filtering \
    ${OUT_DIR}/${PREFIX_M01}_dels.10xCustomFilt.noBND.vcf.gz \
    ${OUT_DIR}/${PREFIX_M01}_large_svs.10xCustomFilt.noBND.vcf.gz \
    ${OUT_DIR}/${PREFIX_M01}_phased_variants.10xCustomFilt.vcf.gz \
    ${OUT_DIR} \
    ${PREFIX_M01} \
    ${REPEAT_ANN} \
    ${HIGH_COPY_BED} \
    ${REF_Ns_BED}
# M20
pass2_filtering \
    ${OUT_DIR}/${PREFIX_M20}_dels.10xCustomFilt.noBND.vcf.gz \
    ${OUT_DIR}/${PREFIX_M20}_large_svs.10xCustomFilt.noBND.vcf.gz \
    ${OUT_DIR}/${PREFIX_M20}_phased_variants.10xCustomFilt.vcf.gz \
    ${OUT_DIR} \
    ${PREFIX_M20} \
    ${REPEAT_ANN} \
    ${HIGH_COPY_BED} \
    ${REF_Ns_BED}
# M29
pass2_filtering \
    ${OUT_DIR}/${PREFIX_M29}_dels.10xCustomFilt.noBND.vcf.gz \
    ${OUT_DIR}/${PREFIX_M29}_large_svs.10xCustomFilt.noBND.vcf.gz \
    ${OUT_DIR}/${PREFIX_M29}_phased_variants.10xCustomFilt.vcf.gz \
    ${OUT_DIR} \
    ${PREFIX_M29} \
    ${REPEAT_ANN} \
    ${HIGH_COPY_BED} \
    ${REF_Ns_BED}

# Merge and concat files at this stage to explore visually in IGV
# DELs
bcftools merge \
    ${OUT_DIR}/${PREFIX_M01}_dels.noRepeatOverlap.noRefNs.vcf.gz \
    ${OUT_DIR}/${PREFIX_M20}_dels.noRepeatOverlap.noRefNs.vcf.gz \
    ${OUT_DIR}/${PREFIX_M29}_dels.noRepeatOverlap.noRefNs.vcf.gz \
    -O z -o ${OUT_DIR}/${PREFIX}_dels_merged.noRepeatOverlap.noRefNs.vcf.gz
# Large SVs
bcftools merge \
    ${OUT_DIR}/${PREFIX_M01}_large_svs.noRepeatOverlap.noRefNs.vcf.gz \
    ${OUT_DIR}/${PREFIX_M20}_large_svs.noRepeatOverlap.noRefNs.vcf.gz \
    ${OUT_DIR}/${PREFIX_M29}_large_svs.noRepeatOverlap.noRefNs.vcf.gz \
    -O z -o ${OUT_DIR}/${PREFIX}_large_svs_merged.noRepeatOverlap.noRefNs.vcf.gz
# Phased variants
bcftools merge \
    ${OUT_DIR}/${PREFIX_M01}_phased_variants.noRepeatOverlap.noRefNs.vcf.gz \
    ${OUT_DIR}/${PREFIX_M20}_phased_variants.noRepeatOverlap.noRefNs.vcf.gz \
    ${OUT_DIR}/${PREFIX_M29}_phased_variants.noRepeatOverlap.noRefNs.vcf.gz \
    -O z -o ${OUT_DIR}/${PREFIX}_phased_variants_merged.noRepeatOverlap.noRefNs.vcf.gz
# Index vcfs
tabix -p vcf ${OUT_DIR}/${PREFIX}_phased_variants_merged.noRepeatOverlap.noRefNs.vcf.gz
tabix -p vcf ${OUT_DIR}/${PREFIX}_dels_merged.noRepeatOverlap.noRefNs.vcf.gz
tabix -p vcf ${OUT_DIR}/${PREFIX}_large_svs_merged.noRepeatOverlap.noRefNs.vcf.gz
# Concatenate files for IGV exploration
bcftools concat --allow-overlaps \
    ${OUT_DIR}/${PREFIX}_phased_variants_merged.noRepeatOverlap.noRefNs.vcf.gz \
    ${OUT_DIR}/${PREFIX}_dels_merged.noRepeatOverlap.noRefNs.vcf.gz \
    ${OUT_DIR}/${PREFIX}_large_svs_merged.noRepeatOverlap.noRefNs.vcf.gz \
    | bcftools sort -O v -o ${OUT_DIR}/${PREFIX}_all_var_filt_concat.vcf

# Exclude sites that differ between 10x Morex and Morex reference
#grep "#" ${OUT_DIR}/${PREFIX}_filtMorexDiffs.vcf > ${OUT_DIR}/${PREFIX}_filtered_no_morex_diffs.vcf
bedtools intersect -v -header -a ${OUT_DIR}/${PREFIX}_filtMorexDiffs.vcf -b ${BED_EXCLUSION_LIST} >> ${OUT_DIR}/${PREFIX}_filtered_no_morex_diffs.vcf

# Exclude non-unique variants (i.e., variants that are present in more than 1 of the mutated lines)
cat ${OUT_DIR}/${PREFIX}_filtered_no_morex_diffs.vcf | vcf-annotate -f ${VCFTOOLS_CUSTOM_FILTER} > ${OUT_DIR}/${PREFIX}_filtered_singleton_ann.vcf
# Create a file containing only singletons
# and remove chrUn variants
grep "#" ${OUT_DIR}/${PREFIX}_filtered_singleton_ann.vcf > ${OUT_DIR}/${PREFIX}_filtered_singletons_only.vcf
grep -v "#" ${OUT_DIR}/${PREFIX}_filtered_singleton_ann.vcf | grep -v "chrUn" | grep "SINGLETON" >> ${OUT_DIR}/${PREFIX}_filtered_singletons_only.vcf

# Pull out homozygous sites only
bcftools view -i 'GT[*]="hom"' ${OUT_DIR}/${PREFIX}_filtered_singletons_only.vcf > ${OUT_DIR}/${PREFIX}_filtered_hom_singletons_only.vcf
