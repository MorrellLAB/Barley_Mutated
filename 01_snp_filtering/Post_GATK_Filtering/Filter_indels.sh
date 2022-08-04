#!/bin/bash

set -e
set -o pipefail

# Dependencies
module load bcftools_PML/1.15.1
module load gatk/4.1.2
module load htslib/1.9
module load python3/3.7.4_anaconda2019.10
module load R/4.0.4
module load bedtools/2.29.2
module unload R/3.4.4-tiff # One of the packages above also load an older version of R that will mess with downstream plotting
export PATH=${PATH}:/panfs/roc/groups/9/morrellp/liux1299/Software

# User provided input arguments
VCF="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Create_HC_Subset/mut8_and_hybrid_barley_raw_variants_indels.vcf"
OUT_PREFIX="mut8_and_hybrid_barley_indels"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered"
# Scratch directory used to store intermediate files that don't need to be kept long term
SCRATCH_DIR="/scratch.global/liux1299/temp_mut8_and_hybrid_barley"

# Thresholds for filtering
# Note: for indels, shouldn't filter on mapping quality related annotations
# see GATK docs: https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
# Filter QD < 2
QD="2"
# Filter FS > 200
FS="200"
# Filter SOR > 3
SOR="3"
# Filter ReadPosRankSum < -20.0
ReadPosRankSum="-20.0"
# Quality cutoff
QUAL_CUTOFF="30"
# Min and Max DP per sample threshold
MIN_DP="5"
MAX_DP="86"
# GQ cutoff per sample
GQ_CUTOFF="6"
# Allele balance filter, minimum deviation +/- threshold
MIN_DEV="0.20"
# Max proportion missing
MAX_MISS="0.30"
# Proportion heterozygous genotypes threshold
HET_PROP="0.1"
# Cutoff for filtering fake polymorphic sites
# These are sites that are actually monomorphic but show up as polymorphic
#   due to missingness and how programs handle that
# We'll filter out sites where the proportion is less than expected for a singleton
#   in the population. The minimum should be a singleton.
# To get this cutoff use: 1 / (number of samples * 2)
# Here, we have 35 total samples in the VCF: 1/70 = 0.01428571 (round down)
FAKE_POLY="0.01"

# List of regions where REF has stretches of N's
REF_Ns_BED="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/stretches_of_Ns/Barley_MorexV3_pseudomolecules_parts_missing.bed"
# Repeat annotations
REPEAT_ANN="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/annotation/HvulgareMorex_702_V3.repeatmasked_assembly_V3.parts.gff3"
# High copy regions (e.g., chloroplasts, mitochondria, rDNA repeats, centromere repeats, etc.)
HIGH_COPY_BED="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/high_copy_regions/Morex_v3_high_copy_uniq.parts.bed"

#-----------------
# Check the out dir and scratch dir exist, if not make them
mkdir -p ${OUT_DIR} ${OUT_DIR}/Intermediates ${SCRATCH_DIR}

function check_filter_stringency() {
    local vcf_file="$1"
    local num_sites="$2"
    if [[ ${num_sites} == 0 ]]; then
        echo "No sites left after filtering. Try using less stringent criteria. File with no sites is: ${vcf_file}. Exiting..." >&2
        exit 8 # If not sites left, error out with message
    fi
}

export -f check_filter_stringency

function count_sites() {
    local vcf="$1"
    local log_file="$2"
    # Get the number of sites left after filtering
    if [[ "${vcf}" == *"gz"* ]]; then
        # We are working with gzipped VCF file
        num_sites=$(zgrep -v "#" ${vcf} | wc -l)
    else
        # We are working with uncompressed vcf
        num_sites=$(grep -v "#" ${vcf} | wc -l)
    fi
    # Append the number of sites remaining to file
    printf "${vcf}\t${num_sites}\n" >> ${log_file}
    check_filter_stringency ${vcf} ${num_sites}
}

export -f count_sites

function RemoveComplexVar() {
    local vcf="$1"
    local accession="$2"
    local out_dir="$3"
    # Keep a record of complex variants we removed
    bcftools filter \
        -i '(STRLEN(REF) >= 2 & STRLEN(ALT) >= 2) | (STRLEN(REF) > 2 | STRLEN(ALT) > 2)' \
        -o ${out_dir}/complex_var_only_${accession}.vcf \
        -O v \
        ${vcf}
    # Create a VCF that excludes the complex variants
    # We want to keep cases like REF is A and ALT is AC, or REF is TC and ALT is T. These are indels.
    # We DO NOT want cases where REF is GTAAA and ALT is G, or REF is A and ALT is ATGC.
    # We also DO NOT want cases where REF is CGCGGGACGAC and ALT is TGCGGGACGAC,T.
    bcftools filter \
        -e '(STRLEN(REF) >= 2 & STRLEN(ALT) >= 2) | (STRLEN(REF) > 2 | STRLEN(ALT) > 2)' \
        -o ${out_dir}/${accession}_noComplex.vcf.gz \
        -O z \
        ${vcf}
    # Index the vcf file
    bcftools index \
        ${out_dir}/${accession}_noComplex.vcf.gz \
        -t \
        -o ${out_dir}/${accession}_noComplex.vcf.gz.tbi
}

export -f RemoveComplexVar

# Get the number of sites left for starting VCF and append to file
count_sites ${VCF} ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

# Remove sites that aren't polymorphic (minor allele count of 0) and unused alternate alleles.
echo "Removing sites that aren't polymorphic and unused alternate alleles..."
gatk SelectVariants \
    -V "${VCF}" \
    --exclude-non-variants true \
    --remove-unused-alternates true \
    -O "${SCRATCH_DIR}/${OUT_PREFIX}_polymorphic.vcf.gz" \
    --tmp-dir ${SCRATCH_DIR}
echo "Done removing sites that aren't polymorphic and unused alternate alleles."
# Get the number of sites left for starting VCF and append to file
count_sites ${SCRATCH_DIR}/${OUT_PREFIX}_polymorphic.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

### Filter on GATK annotations
# Filter out sites with GQ below cutoff, DP between minDP and maxDP, and by allelic balance
# Note: This script sets genotypes to missing if the GQ score is below the cutoff,
#   no sites actually get filtered out at this stage. Sites will get filtered out when we filter on missingness.
#   The same occurs for DP cutoffs, we filter out sites where depth is too low or too high
echo "Removing sites based on the following annotations and values: QD < ${QD}; FS > ${FS}; SOR > ${SOR}, ReadPosRankSum < ${ReadPosRankSum}"
# Exclude sites if any of the following conditions are met
bcftools filter -e "INFO/QD < ${QD} | INFO/FS > ${FS} | INFO/SOR > ${SOR} | INFO/ReadPosRankSum < ${ReadPosRankSum} | QUAL < ${QUAL_CUTOFF}" -O v -o ${SCRATCH_DIR}/${OUT_PREFIX}.filt_QD_FS_SOR_ReadPosRankSum_QUAL.vcf ${SCRATCH_DIR}/${OUT_PREFIX}_polymorphic.vcf.gz
# Get the number of sites left for starting VCF and append to file
count_sites ${SCRATCH_DIR}/${OUT_PREFIX}.filt_QD_FS_SOR_ReadPosRankSum_QUAL.vcf ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

echo "Setting sites to missing based on per sample annotations:"
echo "per sample DP where DP < ${MIN_DP} or DP > ${MAX_DP}"
echo "GQ < ${GQ_CUTOFF}"
# -i in +setGT means if GT ann meet condition, set to missing
bcftools +setGT ${SCRATCH_DIR}/${OUT_PREFIX}.filt_QD_FS_SOR_ReadPosRankSum_QUAL.vcf -- -t q -n "." -i "FMT/DP<${MIN_DP} | FMT/DP > ${MAX_DP} | FMT/GQ<${GQ_CUTOFF}" > ${SCRATCH_DIR}/${OUT_PREFIX}.filt_QD_FS_SOR_ReadPosRankSum_QUAL.DP_GQ.vcf

# Filter out sites with high missingness
echo "Removing sites that missingness exceeding threshold ${MAX_MISS}..."
bcftools filter -e "F_PASS(GT='mis') > ${MAX_MISS}" ${SCRATCH_DIR}/${OUT_PREFIX}.filt_QD_FS_SOR_ReadPosRankSum_QUAL.DP_GQ.vcf -o ${SCRATCH_DIR}/${OUT_PREFIX}.filtMISS.vcf
echo "Done removing sites with missingness exceeding threshold."
# Get the number of sites left after filtering and append to file
count_sites ${SCRATCH_DIR}/${OUT_PREFIX}.filtMISS.vcf ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

# Remove complex variants (>=3 in length) in REF and ALT alleles columns
RemoveComplexVar ${SCRATCH_DIR}/${OUT_PREFIX}.filtMISS.vcf ${OUT_PREFIX} ${OUT_DIR}
# Get the number of sites left after filtering and append to file
count_sites ${OUT_DIR}/${OUT_PREFIX}_noComplex.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

# Remove indels that overlap with repeat annotated regions, high copy regions, and that overlap with stretches of Ns
bedtools intersect -wa -v -header -a ${OUT_DIR}/${OUT_PREFIX}_noComplex.vcf.gz -b ${REPEAT_ANN} ${HIGH_COPY_BED} ${REF_Ns_BED} | bgzip > ${OUT_DIR}/${OUT_PREFIX}_noComplex.noRepeatOverlap.noRefNs.vcf.gz
# Get the number of sites left after filtering and append to file
count_sites ${OUT_DIR}/${OUT_PREFIX}_noComplex.noRepeatOverlap.noRefNs.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_num_sites.log
