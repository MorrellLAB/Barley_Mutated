#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=12gb
#SBATCH --tmp=8gb
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Dependencies
module load bcftools/1.10.2
module load htslib/1.9
module load python3/3.7.4_anaconda2019.10
module load vcftools_ML/0.1.16
module load bedtools/2.29.2
module unload R/3.4.4-tiff # One of the packages above also load an older version of R that will mess with downstream plotting

# User provided input arguments
#VCF="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/mut8_and_hybrid_barley_snps_poly_and_biallelic.vcf.gz"
VCF="/scratch.global/liux1299/temp_mut8_and_hybrid_barley/mut8_and_hybrid_barley_snps_polymorphic.vcf.gz"
OUT_PREFIX="mut8_and_hybrid_barley_snps"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered"
# Scratch directory used to store intermediate files that don't need to be kept long term
SCRATCH_DIR="/scratch.global/liux1299/temp_mut8_and_hybrid_barley"
# Thresholds for filtering
# Proportion heterozygous genotypes threshold
HET_PROP="0.1"
# Min and Max DP per sample threshold
MIN_DP="5"
MAX_DP="158"
# Max proportion missing
MAX_MISS="0.30"
# Quality cutoff
QUAL_CUTOFF="30"
# GQ cutoff per sample
GQ_CUTOFF="9"
# Allele balance filter, minimum deviation +/- threshold
MIN_DEV="0.10"
# Cutoff for filtering fake polymorphic sites
# These are sites that are actually monomorphic but show up as polymorphic
#   due to missingness and how programs handle that
# We'll filter out sites where the proportion is less than expected for a singleton
#   in the population. The minimum should be a singleton.
# To get this cutoff use: 1 / (number of samples * 2)
# Here, we have 35 total samples in the VCF: 1/70 = 0.01428571 (round down)
FAKE_POLY="0.01"

# List of regions where REF has stretches of N's
REF_Ns_BED="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/stretches_of_Ns/Barley_MorexV3_pseudomolecules_parts_missing.bed"
# Repeat annotations
REPEAT_ANN="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/annotation/HvulgareMorex_702_V3.repeatmasked_assembly_V3.parts.gff3"
# High copy regions (e.g., chloroplasts, mitochondria, rDNA repeats, centromere repeats, etc.)
HIGH_COPY_BED="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/high_copy_regions/Morex_v3_high_copy_uniq.parts.bed"

# Define paths to executables and scripts
#FILT_AB_DP_GQ_SCRIPT="/panfs/jay/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/01_snp_filtering/Post_GATK_Filtering/filter_vcf_AB_DP_GQ.py"
FILT_AB_SCRIPT="/panfs/jay/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/01_snp_filtering/Post_GATK_Filtering/filter_vcf_AB.py"

#-----------------
# Check the out dir and scratch dir exist, if not make them
mkdir -p ${OUT_DIR} ${SCRATCH_DIR} ${OUT_DIR}/vcf_summary

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

# Get the number of sites left for starting VCF and append to file
count_sites ${VCF} ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

# Filter out sites with GQ below cutoff, DP between minDP and maxDP, and by allelic balance
# Note: This script sets genotypes to missing if the GQ score is below the cutoff,
#   no sites actually get filtered out at this stage. Sites will get filtered out when we filter on missingness.
#   The same occurs for DP cutoffs, we filter out sites where depth is too low or too high
echo "Removing sites below GQ threshold ${GQ_CUTOFF}..."
echo "Also removing genotypes where per sample DP < ${MIN_DP} or per sample DP > ${MAX_DP}..."
echo "Filtering by allelic balance where deviation is ${MIN_DEV}"
python3 ${FILT_AB_SCRIPT} ${VCF} ${MIN_DEV} > ${SCRATCH_DIR}/${OUT_PREFIX}_AB_filt.vcf
# Use bcftools to set genotypes to missing based on cutoffs
# -i in +setGT means if GT ann meet condition, set to missing
bcftools +setGT ${SCRATCH_DIR}/${OUT_PREFIX}_AB_filt.vcf -- -t q -n "." -i "FMT/DP<${MIN_DP} | FMT/DP > ${MAX_DP} | FMT/GQ<${GQ_CUTOFF}" > ${SCRATCH_DIR}/${OUT_PREFIX}_AB_GQ_DP_filt.vcf
echo "Done setting sites below GQ threshold to missing."
echo "Done removing sites where per sample depth is either too low or too high."
AB_GQ_DP_FILT_VCF="${SCRATCH_DIR}/${OUT_PREFIX}_AB_GQ_DP_filt.vcf"

# Filter out sites with high missingness
echo "Removing sites that missingness exceeding threshold ${MAX_MISS}..."
bcftools filter -e "F_PASS(GT='mis') > ${MAX_MISS}" ${AB_GQ_DP_FILT_VCF} -o ${SCRATCH_DIR}/${OUT_PREFIX}_filtMISS.vcf
MISS_FILT_VCF="${SCRATCH_DIR}/${OUT_PREFIX}_filtMISS.vcf"
echo "Done removing sites with missingness exceeding threshold."
# Get the number of sites left after filtering and append to file
count_sites ${MISS_FILT_VCF} ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

# Filter on heterozygosity (filter out highly heterozygous sites)
echo "Removing sites with more than ${HET_PROP} heterozygous genotypes..."
bcftools filter -i "COUNT(GT='het')/(N_SAMPLES-N_MISSING) < ${HET_PROP}" ${MISS_FILT_VCF} -o ${SCRATCH_DIR}/${OUT_PREFIX}_filtered_het.vcf
echo "Done removing highly heterozygous sites."
HET_FILT_VCF="${SCRATCH_DIR}/${OUT_PREFIX}_filtered_het.vcf"
# Get the number of sites left after filtering and append to file
count_sites ${HET_FILT_VCF} ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

# Filter out sites with quality below cutoff
echo "Removing sites below quality threshold ${QUAL_CUTOFF}..."
bcftools filter -e "QUAL < ${QUAL_CUTOFF}" ${HET_FILT_VCF} -O z -o ${SCRATCH_DIR}/${OUT_PREFIX}_filtQUAL.vcf.gz
echo "Done removing sites below quality thresold."
QUAL_FILT_VCF="${SCRATCH_DIR}/${OUT_PREFIX}_filtQUAL.vcf.gz"
# Index vcf
tabix -p vcf ${QUAL_FILT_VCF}
# Get the number of sites left after filtering and append to file
count_sites ${QUAL_FILT_VCF} ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

# Filter out fake polymorphic sites
# These are sites that are actually monomorphic but got labeled as polymorphic
# This is just an extra check to be safe
echo "Removing fake polymorphic sites below MAF threshold ${FAKE_POLY}..."
#vcftools --gzvcf ${QUAL_FILT_VCF} --maf ${FAKE_POLY} --recode --recode-INFO-all --out ${OUT_DIR}/${OUT_PREFIX}_final
vcftools --gzvcf ${QUAL_FILT_VCF} --maf ${FAKE_POLY} --recode --recode-INFO-all --out ${OUT_DIR}/${OUT_PREFIX}_noFakePoly
echo "Done removing fake polymorphic sites."
# Rename output and compress
mv ${OUT_DIR}/${OUT_PREFIX}_noFakePoly.recode.vcf ${OUT_DIR}/${OUT_PREFIX}_noFakePoly.vcf
bgzip ${OUT_DIR}/${OUT_PREFIX}_noFakePoly.vcf
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}_noFakePoly.vcf.gz

# Filter to only biallelic sites
echo "Filter to only biallelic sites..."
bcftools view -m2 -M2 ${OUT_DIR}/${OUT_PREFIX}_noFakePoly.vcf.gz -O z -o "${OUT_DIR}/${OUT_PREFIX}_biallelic.vcf.gz"
# Index VCF
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}_biallelic.vcf.gz
# Get the number of sites left after filtering and append to file
count_sites ${OUT_DIR}/${OUT_PREFIX}_biallelic.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

# Pull out triallelic sites (for exploratory purposes only)
echo "Pulling out triallelic sites..."
bcftools view -m3 ${OUT_DIR}/${OUT_PREFIX}_noFakePoly.vcf.gz -O z -o "${OUT_DIR}/${OUT_PREFIX}_triallelic.vcf.gz"

# Remove variants that overlap with repeat annotated regions, high copy regions, and that overlap with stretches of Ns
bedtools intersect -wa -v -header -a ${OUT_DIR}/${OUT_PREFIX}_biallelic.vcf.gz -b ${REPEAT_ANN} ${HIGH_COPY_BED} ${REF_Ns_BED} | bgzip > ${OUT_DIR}/${OUT_PREFIX}_biallelic.noRepeatOverlap.noRefNs.vcf.gz
# Get the number of sites left after filtering and append to file
count_sites ${OUT_DIR}/${OUT_PREFIX}_biallelic.noRepeatOverlap.noRefNs.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_num_sites.log
