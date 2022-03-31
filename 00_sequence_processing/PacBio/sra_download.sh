#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=22gb
#SBATCH --tmp=12gb
#SBATCH -t 05:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Dependencies
module load parallel/20210822

# User provided input arguments
FTP_LIST="/panfs/roc/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/00_sequence_processing/PacBio/sra_ftp_links.txt"
OUT_DIR="/scratch.global/liux1299/sra_pacbio"

#---------------------
# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR}
# Go into OUT_DIR
cd ${OUT_DIR}

# Download sra files
parallel --verbose wget {} :::: ${FTP_LIST}
