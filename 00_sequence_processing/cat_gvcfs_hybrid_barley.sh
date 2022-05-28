#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=12gb
#SBATCH --tmp=8gb
#SBATCH -t 02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

set -e
set -o pipefail

# This script concatenates GVCFs split into regions for each sample and utilizes
#   Slurm job arrays for parallelization across samples.

# Dependencies
module load java/openjdk-8_202
PICARD_JAR=/panfs/roc/groups/9/morrellp/public/Software/picard_ML_2.23.1/picard.jar

# User provided input arguments
# List of GVCF files
GVCF_LIST="/scratch.global/liux1299/hybrid_barley_morex_v3/Haplotype_Caller/gvcf_list.txt"
SAMPLE_NAMES="/scratch.global/liux1299/hybrid_barley_morex_v3/Haplotype_Caller/hybrid_barley_sample_names.txt"
OUT_DIR="/scratch.global/liux1299/hybrid_barley_morex_v3/gvcf_concat"
TEMP="/scratch.global/liux1299/hybrid_barley_morex_v3/gvcf_concat/temp"
MEM="10g"

#--------------
mkdir -p ${OUT_DIR} ${TEMP}

# Prepare array for utilizing Slurm job arrays
SAMP_ARR=($(cat ${SAMPLE_NAMES}))
# Determine maximum array limit
MAX_ARRAY_LIMIT=$[${#SAMP_ARR[@]} - 1]
echo "Maximum array limit is ${MAX_ARRAY_LIMIT}."

# Get the current sample we are processing
CURR_SAMPLE=${SAMP_ARR[${SLURM_ARRAY_TASK_ID}]}
echo "Currently processing sample: ${CURR_SAMPLE}"

# Prepare list of GVCF regions for current sample we are concatenating
# Input list to Picard must end in file extension .list
grep "${CURR_SAMPLE}" "${GVCF_LIST}" | sort -V > ${OUT_DIR}/chr_split_gvcf_list-${CURR_SAMPLE}.list

# Concatenate split regions for current sample into single GVCF
java -Xmx"${MEM}" -jar ${PICARD_JAR} GatherVcfs \
    INPUT="${OUT_DIR}/chr_split_gvcf_list-${CURR_SAMPLE}.list" \
    OUTPUT="${OUT_DIR}/${CURR_SAMPLE}_RawGLs.g.vcf" \
    VERBOSITY="WARNING" \
    TMP_DIR="${TEMP}"
