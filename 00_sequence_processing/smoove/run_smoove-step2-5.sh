#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=56gb
#SBATCH --tmp=40gb
#SBATCH -t 36:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# This script runs smoove for large cohorts (n >= ~ 40)
#   Note: cohort size in smoove seems to be based on the human genome size

# Dependencies
module load parallel/20210822
module load samtools/1.9
module load htslib/1.9
module load bcftools/1.10.2
# svtools and svtyper both require python 2.7
module load python2/2.7.16_anaconda2019.10
source activate /panfs/roc/groups/9/morrellp/liux1299/.conda/envs/my_svtyper
# Export path to smoove
export PATH=${PATH}:/panfs/roc/groups/9/morrellp/shared/Software/smoove

# User provided input arguments
# Directory containing *smoove.genotyped.vcf.gz files for each sample
GENO_VCF_DIR="/scratch.global/liux1299/results_smoove/called_genotypes"
# List of BAM files
BAM_LIST="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/smoove_processing/bam_symlinks/wgs_bam_list.txt"
OUT_NAME="mut_barley_cohort"
REF_FASTA="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta"
# Full filepath to GFF file
GFF="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/gene_annotation/Hv_Morex.pgsb.Jul2020.sorted.parts.gff3"
# Output final output files here
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/smoove_processing/results_smoove"
# Output intermediate files to temporary storage
SCRATCH_DIR="/scratch.global/liux1299/results_smoove"

#----------------------
# Make output directories
mkdir -p ${OUT_DIR} ${OUT_DIR}/results_genotyped ${SCRATCH_DIR} ${SCRATCH_DIR}/temp ${SCRATCH_DIR}/merged_sites ${SCRATCH_DIR}/results_genotyped_sample
# For large cohorts, set temp dir
export TMPDIR="${SCRATCH_DIR}/temp"

# Run smoove
# Step 2: Get the union of sites across all samples (can parallelize this across as many CPUs or machines as needed):
# this will create ${SCRATCH_DIR}/merged_sites/merged.sites.vcf.gz
echo "Merging sites..."
smoove merge \
    --name ${OUT_NAME}_merged \
    --fasta ${REF_FASTA} \
    --outdir ${SCRATCH_DIR}/merged_sites \
    ${GENO_VCF_DIR}/*.genotyped.vcf.gz
echo "Done merging sites."

# Step 3: Genotype each sample at those sites (this can parallelize this across as many CPUs or machines as needed) and run duphold to add depth annotations.
function smoove_genotype() {
    local curr_bam="$1"
    local vcf_merged="$2"
    local out_dir="$3"
    local ref_fasta="$4"
    # Check if sample name contains substring
    if [[ "${curr_bam}" == *"_phased_possorted_bam.bam"* ]]
    then
        # This is a 10x Genomics sample, generate clean sample name accordingly
        sample_name=$(basename ${curr_bam} _phased_possorted_bam.bam)
    else
        sample_name=$(basename ${curr_bam} .bam)
    fi
    # Genotype each sample at each one of those sites to extract population-aware genotyping calls
    smoove genotype \
        --duphold \
        --removepr \
        --processes 2 \
        --name ${sample_name} \
        --outdir ${out_dir} \
        --fasta ${ref_fasta} \
        --vcf ${vcf_merged} \
        ${curr_bam}
}

export -f smoove_genotype

echo "Genotyping samples..."
parallel --verbose smoove_genotype {} ${SCRATCH_DIR}/merged_sites/${OUT_NAME}_merged.sites.vcf.gz ${SCRATCH_DIR}/results_genotyped_sample ${REF_FASTA} :::: ${BAM_LIST}
echo "Done genotyping samples."

# Step 4: Paste all the single sample VCFs with the same number of variants to get a single, squared, joint-called file.
#   The VCFs should be a list of VCFs output from smoove genotype
echo "Pasting single sample VCFs..."
smoove paste \
    --name ${OUT_NAME}_sites \
    --outdir ${OUT_DIR}/results_genotyped \
    ${SCRATCH_DIR}/results_genotyped_sample/*.vcf.gz
echo "Done pasting single sample VCFs."

# Step 5: (optional) Annotate the variants with exons, UTRs that overlap from a GFF and annotate high-quality heterozygotes:
echo "Annotating variants from GFF..."
smoove annotate --gff ${GFF} ${OUT_DIR}/results_genotyped/${OUT_NAME}_sites.smoove.square.vcf.gz | bgzip -c > ${OUT_DIR}/results_genotyped/${OUT_NAME}_sites.smoove.square.anno.vcf.gz
echo "Done annotating variants from GFF."
