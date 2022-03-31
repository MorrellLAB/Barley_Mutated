#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=12gb
#SBATCH --tmp=8gb
#SBATCH -t 02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# This calls SVs using Sniffles2

# Dependencies
module load samtools/1.9

# User provided input arguments
BAM="/panfs/roc/groups/9/morrellp/liux1299/My_IGV/Morex_1_14_align_V2_sorted_parts.bam"

#-------------------
# Index BAM
samtools index -b ${BAM}
