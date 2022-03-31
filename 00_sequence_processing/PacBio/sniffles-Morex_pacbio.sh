#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=36gb
#SBATCH --tmp=22gb
#SBATCH -t 05:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Call SVs using Sniffles2

# Dependencies
module load samtools/1.9
module load python3/3.8.3_anaconda2020.07_mamba
# Activate conda environment
source activate /panfs/roc/groups/9/morrellp/liux1299/.conda/envs/sniffles_env

# User provided input arguments
BAM="/scratch.global/liux1299/sra_pacbio/pacbio_morex_v2/Morex_pacbio/Morex_pacbio_90.bam"
REF="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules_parts.fasta"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/pacbio_morex/pacbio_morex_v2"
NUM_THREADS="16"

#-----------------
# Make output directory
mkdir -p ${OUT_DIR}

# Check if BAM file is indexed, if not index it
# Note: this may need to be changed to CSI indexing depending on
#   if we are dealing with the parts reference or pseudomolecular reference
if ! [ -f "${BAM}.bai" ]; then
    # Generate BAI index
    echo "Generating BAI index for BAM file..."
    samtools index -b ${BAM}
fi

OUT_PREFIX=$(basename ${BAM} .bam)

# Output VCF
sniffles --input ${BAM} \
    --vcf ${OUT_DIR}/${OUT_PREFIX}.vcf.gz \
    --reference ${REF} \
    --threads ${NUM_THREADS}
