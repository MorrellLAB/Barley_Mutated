#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=12gb
#SBATCH --tmp=8gb
#SBATCH -t 04:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# This script removes phase info from GATK VCF files for downstream tools where
#   phasing causes an error in the software run.

# Dependencies
module load python3/3.8.3_anaconda2020.07_mamba
source activate /panfs/jay/groups/9/morrellp/liux1299/.conda/envs/whatshap_env

# User provided input arguments
VCF1="/scratch.global/liux1299/temp_mut8_and_hybrid_barley/mut8_and_hybrid_barley_snps_polymorphic.vcf.gz"
# After filtering on various annotations
VCF2="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/mut8_and_hybrid_barley_snps_biallelic.vcf.gz"
VCF3="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/mut8_and_hybrid_barley_snps_biallelic.noRepeatOverlap.noRefNs.vcf.gz"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered/unphased_vcf"

#----------------------
mkdir -p ${OUT_DIR}

function unphase_genotypes() {
    local vcf="$1"
    local out_dir="$2"
    # Get out_prefix
    if [[ "${vcf}" == *"gz"* ]]; then
        # We are working with gzipped VCF file
        out_prefix=$(basename ${vcf} .vcf.gz)
    else
        # We are working with uncompressed vcf
        out_prefix=$(basename ${vcf} .vcf)
    fi
    whatshap unphase ${vcf} > ${out_dir}/${out_prefix}_unphased.vcf
}

export -f unphase_genotypes

unphase_genotypes ${VCF3} ${OUT_DIR}
unphase_genotypes ${VCF2} ${OUT_DIR}
unphase_genotypes ${VCF1} ${OUT_DIR}
