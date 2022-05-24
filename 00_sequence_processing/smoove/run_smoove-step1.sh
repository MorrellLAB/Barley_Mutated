#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=56gb
#SBATCH --tmp=40gb
#SBATCH -t 40:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

set -e
set -o pipefail

# This script runs smoove for large cohorts (n >= ~ 40) STEP 1
#   Note: cohort size in smoove seems to be based on the human genome size
#   STEP 1 can be parallelize across samples, so we'll submit this part as a job array
#       one array index per sample

# Dependencies
module load samtools/1.9
module load htslib/1.9
module load bcftools/1.10.2
# svtools and svtyper both require python 2.7
module load python2/2.7.16_anaconda2019.10
source activate /panfs/roc/groups/9/morrellp/liux1299/.conda/envs/my_svtyper
# Export path to smoove
export PATH=${PATH}:/panfs/roc/groups/9/morrellp/shared/Software/smoove

# User provided input arguments
# List of BAM files
BAM_LIST="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/smoove_processing/bam_symlinks/wgs_bam_list.txt"
REF_FASTA="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta"
# Per-sample exclusion regions file
# Each sample bed includes sample specific high cov regions and
#   regions where there are stretches of N's in the reference genome
EXCLUDE_BED_LIST="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/smoove_processing/masked_regions/per_sample_gap_and_high_depth/per_sample_exclude_bed_list.txt"
# Exclude chrUn (this is formatted as a regex expression according to smoove's documentation)
# For more than one, provide a comma delimited list of exclusions
#   See "smoove call -h" the "--excludechroms" description for formatting help
# Here, an alternative way to format the input is to create a file with exclusion chromosomes,
#   one per line and use the following code to format:
#   EXCLUDE_CHR=$(tr '\n' ',' < exclude_chr_list.txt | sed 's/,$//g')
EXCLUDE_CHR="~^chrUn"
# Output intermediate files to temporary storage
SCRATCH_DIR="/scratch.global/liux1299/results_smoove"

#----------------------
# Make output directories
mkdir -p ${SCRATCH_DIR} ${SCRATCH_DIR}/temp ${SCRATCH_DIR}/called_genotypes
# For large cohorts, set temp dir
export TMPDIR="${SCRATCH_DIR}/temp"

# Prepare array for utilizing Slurm job arrays
BAM_ARR=($(cat ${BAM_LIST}))
# Determine maximum array limit
MAX_ARRAY_LIMIT=$[${#BAM_ARR[@]} - 1]
echo "Maximum array limit is ${MAX_ARRAY_LIMIT}."

# Get the current BAM file we are processing
CURR_BAM=${BAM_ARR[${SLURM_ARRAY_TASK_ID}]}
echo "Currently processing bam file: ${CURR_BAM}"

# Define functions for parallelization
# Function to run Step 1
function smoove_call_genotypes() {
    local curr_sample="$1"
    local out_dir="$2"
    local exclude_bed_list="$3"
    local ref_fasta="$4"
    local exclude_chr="$5"
    # Check if sample name contains substring
    if [[ "${curr_sample}" == *"_phased_possorted_bam.bam"* ]]
    then
        # This is a 10x Genomics sample, generate clean sample name accordingly
        sample_name=$(basename ${curr_sample} _phased_possorted_bam.bam)
    else
        sample_name=$(basename ${curr_sample} .bam)
    fi
    # Get sample specific exclusion regions bed file
    samp_exclude_regions=$(grep ${sample_name} ${exclude_bed_list})
    # Run smoove call
    smoove call \
        --outdir ${out_dir} \
        --exclude ${samp_exclude_regions} \
        --excludechroms ${exclude_chr} \
        --name ${sample_name} \
        --fasta ${ref_fasta} \
        --processes 2 \
        --genotype ${curr_sample}
}

export -f smoove_call_genotypes

# Run smoove
# Step 1: For each sample, call genotypes
smoove_call_genotypes ${CURR_BAM} ${SCRATCH_DIR}/called_genotypes ${EXCLUDE_BED_LIST} ${REF_FASTA} ${EXCLUDE_CHR}
