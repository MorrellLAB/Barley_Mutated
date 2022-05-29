#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=6gb
#SBATCH --tmp=4gb
#SBATCH -t 00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

set -e
set -o pipefail

module load bcftools/1.9

# Removes complex variants (>=3 in length) in REF and ALT alleles columns

# User provided arguments
# List of full filepaths to phased variants VCF files
VCF_LIST="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/phased_var_vcf_list.txt"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/complex_variants"

#----------------------
mkdir -p ${OUT_DIR}

# Prepare array for utilizing Slurm job arrays
VCF_ARR=($(cat ${VCF_LIST}))
# Determine maximum array limit
MAX_ARRAY_LIMIT=$[${#VCF_ARR[@]} - 1]
echo "Maximum array limit is ${MAX_ARRAY_LIMIT}."

# Get the current VCF file we are processing
CURR_VCF=${VCF_ARR[${SLURM_ARRAY_TASK_ID}]}
echo "Currently processing VCF file: ${CURR_VCF}"
# Get the current sample name
SAMPLE_NAME=$(basename ${CURR_VCF} _phased_variants.vcf.gz)

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
        -o ${out_dir}/${accession}_phased_variants_noComplex.vcf.gz \
        -O z \
        ${vcf}
    # Index the vcf file
    bcftools index \
        ${out_dir}/${accession}_phased_variants_noComplex.vcf.gz \
        -t \
        -o ${out_dir}/${accession}_phased_variants_noComplex.vcf.gz.tbi
}

export -f RemoveComplexVar

# Run function on each VCF file
RemoveComplexVar ${CURR_VCF} ${SAMPLE_NAME} ${OUT_DIR}
