#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=22gb
#SBATCH --tmp=18gb
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Dependencies
module load longranger/2.2.2

# User provided input arguments
# Full filepath to the reference genome fasta file
REF_FASTA=

#-----------------------
# Create compatible reference reference
longranger mkref ${REF_FASTA}
