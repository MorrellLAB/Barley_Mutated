#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=12gb
#SBATCH --tmp=10gb
#SBATCH -t 10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

set -e
set -o pipefail

# This script utilizes Slurm job arrays to run GATK Genotype GVCFs to output invariant sites

# Dependencies
module load java/openjdk-8_202
module load gatk/4.1.2

# User provided input arguments
GVCF_LIST="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/Haplotype_Caller/gvcf_list.txt"
# This should be the same list as used to parallelize across regions during Haplotype Caller
CHR_LIST="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/chr_parts_list.intervals"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/Genotype_GVCFs"
REF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta"

#   What is the nucleotide diversity per base pair (Watterson's theta)?
#       For barley: 0.008
#       For soybean: 0.001
THETA="0.008"

memory="9g"

# Temporary directory
TEMP="/scratch.global/liux1299/gatk_temp"

#--------------
mkdir -p ${OUT_DIR} ${TEMP}

chr_arr=($(cat ${CHR_LIST}))
# Slurm job array
curr_chr="${chr_arr[${SLURM_ARRAY_TASK_ID}]}"
echo "Current chromosome we are processing is: ${curr_chr}"

current_gvcf=$(grep "${curr_chr}" ${GVCF_LIST})
echo "Current GVCF file we are processing is: ${current_gvcf}"

# Genotype GVCFs for single sample
gatk --java-options "-Xmx${memory}" GenotypeGVCFs \
    --reference "${REF}" \
    --variant ${current_gvcf} \
    --include-non-variant-sites \
    -L ${curr_chr} \
    --heterozygosity "${THETA}" \
    --output "${OUT_DIR}/${curr_chr}.vcf.gz" \
    --tmp-dir ${TEMP}
