#!/bin/bash

set -e
set -o pipefail

# Evaluate filtered 10x Genomics VCF file containing 3 mutated lines (includes SNPs, dels, and SVs)

# Dependencies
module load gatk/4.1.2
module load bcftools/1.9
# Note: Picard Jar filepath will need to be modified if not running on MSI
PICARD_JAR=/home/morrellp/public/Software/picard_ML_2.20.2/picard.jar

# User provided input arguments
# Using pass2 vcf because in the script variant_filtering-mut_3_lines.sh, this is
# the file that gets created after filtering based on quality, depth, etc.
VCF=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated/Filtered/mut_3_lines_filtered_pass2.vcf
# Known variants VCF file
KNOWN_VCF=/panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/50k_morex_v2_idt90_parts.vcf
# Reference dict
REF_DICT=/home/morrellp/liux1299/Shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules_parts.dict
# Output file prefix
OUT_PREFIX=mut_3_lines_50k
# Full path to output dir
OUT_DIR=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated/Filtered/Evaluate_Filtered_Callset

# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR}

function index_vcf() {
    local vcf="$1"
    # Check if VCF is indexed, if not index it
    if [ -n "$(ls -A ${vcf}.idx 2>/dev/null)" ]; then
        echo "VCF is already indexed, proceed to next step."
    else
        echo "Indexing ${vcf} file..."
        gatk IndexFeatureFile -F ${vcf}
        echo "Finished indexing ${vcf} file."
    fi
}

export -f index_vcf

# Fix VCF header contigs so that they have a length field
VCF_PREFIX=$(basename ${VCF} .vcf)
java -jar ${PICARD_JAR} UpdateVcfSequenceDictionary \
    I=${VCF} \
    O=${OUT_DIR}/${VCF_PREFIX}_fixedHeader.vcf \
    SD=${REF_DICT}

# Check if vcf files are indexed, if not index it
index_vcf ${OUT_DIR}/${VCF_PREFIX}_fixedHeader.vcf
index_vcf ${KNOWN_VCF}

# The 10x Genomics VCF format allows the end position of SVs to exceed the contig length
# and this causes issues when using GATK tools to evaluate the filtering. So, I will temporarily
# remove these variants for the purposes of evaluating filtering criteria.
#   These problem variants were identified by running:
#   gatk ValidateVariants -R /home/morrellp/liux1299/Shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules_parts.fasta -V mut_3_lines_filtered_pass2_fixedHeader.vcf
# Prepare header lines
grep "#" ${OUT_DIR}/mut_3_lines_filtered_pass2_fixedHeader.vcf > ${OUT_DIR}/temp_excluded_from_eval.vcf
# Identify all positions where the END (stop) exceeds the contig length
# Example: chr1H_part1 start:205170000 stop:205510000
# Create array storing chromosome and chr length
chr_arr=($(grep "##contig" ${OUT_DIR}/temp_excluded_from_eval.vcf | sed -e 's/##contig=<ID=//' -e 's/>//' -e 's/length=//'))

awk '$2=="205170000"' ${OUT_DIR}/mut_3_lines_filtered_pass2_fixedHeader.vcf >> ${OUT_DIR}/temp_excluded_from_eval.vcf
bcftools view -i'END<"${test}"' temp_excluded_from_eval.vcf

# Remove exclusion sites temporarily
bcftools filter -T ^${OUT_DIR}/temp_excluded_from_eval.vcf \
    -o ${OUT_DIR}/mut_3_lines_filtered_pass2_valid_sites.vcf \
    -O v \
    ${OUT_DIR}/mut_3_lines_filtered_pass2_fixedHeader.vcf
# Index VCF
index_vcf ${OUT_DIR}/mut_3_lines_filtered_pass2_valid_sites.vcf

# Compare callset against known population callset
gatk CollectVariantCallingMetrics \
     -I ${OUT_DIR}/mut_3_lines_filtered_pass2_valid_sites.vcf \
     --DBSNP ${KNOWN_VCF} \
     -SD ${REF_DICT} \
     -O ${OUT_DIR}/${OUT_PREFIX}_metrics
