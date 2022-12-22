#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=22gb
#SBATCH --tmp=12gb
#SBATCH -t 10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

set -e
set -o pipefail

# Concatenate variant and invariant sites
# Run as Slurm job array because invariant sites file is too large
# We'll have to concatenate and sort by chromosome parts first, then concatenate
#   sorted chromosome parts together
# Array index will depend on the number of chromosomes (or chromosome parts)

# Dependencies
module load bcftools/1.10.2
module load htslib/1.9

# User provided input arguments
# Variant sites
VCF_VARIANT="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/pixy_input/morex-sample2_phased_variants-snps.DPfilt.noRepeatOverlap.noRefNs.concordant.vcf.gz"
# Invariant sites
VCF_INVARIANT="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/pixy_input/morex-sample2.invariant_sites.filtDP-RGQ-QUAL.snps.concordant.vcf.gz"
# List of chromosome part names
CHR_LIST="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/chr_parts_list.intervals"
# Output file prefix
OUT_PREFIX="morex-sample2_snps_filt_concordant.all_sites"

# Temporary output directory
TEMP_DIR="/scratch.global/liux1299/temp_morex-sample2"

MEMORY="12G"

# Temporary directory for sort
export TMPDIR="/scratch.global/liux1299"

#-------------------
function concat_and_sort() {
    local vcf_variant="$1"
    local vcf_invariant="$2"
    local curr_chr="$3"
    local temp_dir="$4"
    local out_prefix="$5"
    #local memory="$6"
    # Prepare variant vcf output file prefix
    if [[ "${vcf_variant}" == *"gz"* ]]; then
        # We are working with gzipped VCF file
        bn_variant=$(basename ${vcf_variant} .vcf.gz)
    else
        # We are working with uncompressed vcf
        bn_variant=$(basename ${vcf_variant} .vcf)
    fi
    # Prepare invariant vcf output file prefix
    if [[ "${vcf_invariant}" == *"gz"* ]]; then
        # We are working with gzipped VCF file
        bn_invariant=$(basename ${vcf_invariant} .vcf.gz)
    else
        # We are working with uncompressed vcf
        bn_invariant=$(basename ${vcf_invariant} .vcf)
    fi
    # Select current chromosome
    # Variant sites
    bcftools view --regions "${curr_chr}" ${vcf_variant} -O v -o ${temp_dir}/${curr_chr}-${bn_variant}.vcf
    # Invariant sites
    bcftools view --regions "${curr_chr}" ${vcf_invariant} -O v -o ${temp_dir}/${curr_chr}-${bn_invariant}.vcf
    # Concatenate and sort current chromosome sites
    #bcftools concat ${temp_dir}/${curr_chr}-${bn_invariant}.vcf ${temp_dir}/${curr_chr}-${bn_variant}.vcf | bcftools sort --temp-dir ${temp_dir} --max-mem ${memory} -O v -o ${temp_dir}/${curr_chr}-${out_prefix}.vcf
    bcftools concat ${temp_dir}/${curr_chr}-${bn_invariant}.vcf ${temp_dir}/${curr_chr}-${bn_variant}.vcf -O v -o ${temp_dir}/temp_${curr_chr}-${out_prefix}.vcf
    # Sort VCF
    # Use Unix sort since we kept having errors with bcftools sort
    grep "#" ${temp_dir}/temp_${curr_chr}-${out_prefix}.vcf > ${temp_dir}/${curr_chr}-${out_prefix}.vcf
    grep -v "#" ${temp_dir}/temp_${curr_chr}-${out_prefix}.vcf | sort -k1,1 -k2,2n >> ${temp_dir}/${curr_chr}-${out_prefix}.vcf
    # Remove intermediate file
    rm ${temp_dir}/temp_${curr_chr}-${out_prefix}.vcf
}

export -f concat_and_sort

mkdir -p ${OUT_DIR} ${TEMP_DIR}

chr_arr=($(cat ${CHR_LIST}))
# Prepare Slurm job array index
curr_chr="${chr_arr[${SLURM_ARRAY_TASK_ID}]}"
echo "Current chromosome we are processing is: ${curr_chr}"

concat_and_sort ${VCF_VARIANT} ${VCF_INVARIANT} ${curr_chr} ${TEMP_DIR} ${OUT_PREFIX}
