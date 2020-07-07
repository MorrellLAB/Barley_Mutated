#!/bin/bash
#PBS -l mem=8gb,nodes=1:ppn=1,walltime=24:00:00
#PBS -m abe
#PBS -M liux1299@umn.edu
#PBS -q mesabi

set -e
set -o pipefail

# Dependencies
module load python3/3.6.3_anaconda5.0.1
module load parallel/20180922
# Directory where custom script is located
export PATH=${PATH}:/panfs/roc/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/01_snp_filtering/vcf_exploration

# User provided arguments
VCF_LIST=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated/mut_vcf_list.txt

function vcf2h5() {
    local vcf=$1
    # Convert VCF to hdf5
    vcf_to_h5.py ${vcf}
}

export -f vcf2h5

# Store vcf list in array
VCF_ARR=($(cat ${VCF_LIST}))

# Run program
parallel vcf2h5 {} ::: ${VCF_ARR[@]}
