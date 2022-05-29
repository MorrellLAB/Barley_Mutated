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

# Removes mismatching REF allele positions using Slurm job arrays.

# User provided arguments
# IMPORTANT: Make sure vcf list and log list are sorted in the same order
# List of large_svs vcfs
VCF_LIST="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/largeSVs_vcf_list.txt"
# List of log files containing REF_MISMATCH identifiers
# See CL README on Github for how log files were created, these files should be named something
# like: temp_ref_check_warn_M01-3-3_large_svs.log
LOG_LIST="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/largeSVs_ref_mismatch_targets_list.txt"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/fix_ref_mismatch"

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
# Get the current REF_MISMATCH targets from the log file list
SAMPLE_NAME=$(basename ${CURR_VCF} _large_svs.vcf.gz)
TARGETS=$(grep ${SAMPLE_NAME} ${LOG_LIST})

# Create targets file of REF_MISMATCH from log files
grep "REF_MISMATCH" ${TARGETS} | cut -f 2,3 > ${OUT_DIR}/ref_mismatch_${SAMPLE_NAME}_large_svs_targets.txt

# Remove positions where REF allele mismatches
bcftools filter \
    -T ^${OUT_DIR}/ref_mismatch_${SAMPLE_NAME}_large_svs_targets.txt \
    -o ${OUT_DIR}/${SAMPLE_NAME}_large_svs_noRefMismatch.vcf.gz \
    -O z \
    ${CURR_VCF}

# Index vcf files
bcftools index \
    ${OUT_DIR}/${SAMPLE_NAME}_large_svs_noRefMismatch.vcf.gz \
    -t \
    -o ${OUT_DIR}/${SAMPLE_NAME}_large_svs_noRefMismatch.vcf.gz.tbi

# Now, pull out only positions where there is a REF allele mismatch (for exploration purposes)
bcftools filter \
    -T ${OUT_DIR}/ref_mismatch_${SAMPLE_NAME}_large_svs_targets.txt \
    -o ${OUT_DIR}/${SAMPLE_NAME}_large_svs_RefMismatchOnly.vcf \
    -O v \
    ${CURR_VCF}
