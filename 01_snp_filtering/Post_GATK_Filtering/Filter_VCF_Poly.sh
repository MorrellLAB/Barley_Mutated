#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=9gb
#SBATCH --tmp=8gb
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Dependencies
module load bcftools/1.10.2
module load gatk/4.1.2
module load htslib/1.9

# User provided input arguments
VCF="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Variant_Recalibrator/mut8_and_hybrid_barley_snps.recalibrated.pass_sites.vcf.gz"
OUT_PREFIX="mut8_and_hybrid_barley_snps"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/mut8_and_hybrid_barley/Filtered"
# Scratch directory used to store intermediate files that don't need to be kept long term
SCRATCH_DIR="/scratch.global/liux1299/temp_mut8_and_hybrid_barley"

#-----------------
# Check the out dir and scratch dir exist, if not make them
mkdir -p ${OUT_DIR} ${SCRATCH_DIR}

# Check the number of sites in starting VCF
num_sites_vcf=$(zgrep -v "#" ${VCF} | wc -l)
# Append the number of sites remaining to file
printf "${VCF}\t${num_sites_vcf}\n" >> ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

# Remove sites that aren't polymorphic (minor allele count of 0) and unused alternate alleles.
echo "Removing sites that aren't polymorphic and unused alternate alleles..."
gatk SelectVariants \
    -V "${VCF}" \
    --exclude-non-variants true \
    --remove-unused-alternates true \
    -O "${SCRATCH_DIR}/${OUT_PREFIX}_polymorphic.vcf.gz" \
    --tmp-dir ${SCRATCH_DIR}
echo "Done removing sites that aren't polymorphic and unused alternate alleles."
# Get the number of sites left after filtering
num_sites_poly=$(zgrep -v "#" ${SCRATCH_DIR}/${OUT_PREFIX}_polymorphic.vcf.gz | wc -l)
# Append the number of sites remaining to file
printf "${SCRATCH_DIR}/${OUT_PREFIX}_polymorphic.vcf.gz\t${num_sites_poly}\n" >> ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

# # Filter to only biallelic sites
# echo "Filter to only biallelic sites..."
# bcftools view -m2 -M2 "${SCRATCH_DIR}/${OUT_PREFIX}_polymorphic.vcf.gz" -O z -o "${OUT_DIR}/${OUT_PREFIX}_poly_and_biallelic.vcf.gz"
# # Get the number of sites left after filtering
# num_sites_biallelic=$(zgrep -v "#" ${OUT_DIR}/${OUT_PREFIX}_poly_and_biallelic.vcf.gz | wc -l)
# # Append the number of sites remaining to file
# printf "${OUT_DIR}/${OUT_PREFIX}_poly_and_biallelic.vcf.gz\t${num_sites_biallelic}\n" >> ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

# # Index VCF
# tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}_poly_and_biallelic.vcf.gz

# # Pull out triallelic sites (for exploratory purposes only)
# echo "Pulling out triallelic sites..."
# bcftools view -m3 "${SCRATCH_DIR}/${OUT_PREFIX}_polymorphic.vcf.gz" -O z -o "${OUT_DIR}/${OUT_PREFIX}_poly_and_triallelic.vcf.gz"
