#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=90gb
#SBATCH --tmp=70gb
#SBATCH -t 08:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p ram256g,ram1t,amd2tb
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Dependencies
module load python3/3.8.3_anaconda2020.07_mamba
export PATH=${PATH}:/panfs/roc/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/01_snp_filtering/Post_GATK_Filtering

VCF="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Create_HC_Subset/mut8_and_hybrid_barley_raw_variants_indels.vcf"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Create_HC_Subset/Plots"
OUT_PREFIX="mut8_and_hybrid_barley_raw_variants_indels"

#--------------------
# Plot annotations
VCF_ann_viz-indels.py ${VCF} ${OUT_DIR} ${OUT_PREFIX}
