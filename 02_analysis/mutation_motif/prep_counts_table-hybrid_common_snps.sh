#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=10gb
#SBATCH --tmp=8gb
#SBATCH -t 00:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Dependencies
module load bcftools/1.10.2
module load parallel/20210822
module load python3/3.8.3_anaconda2020.07_mamba
module load bedtools/2.29.2
# Conda environment for mutation motif software
source activate /panfs/jay/groups/9/morrellp/liux1299/.conda/envs/mut_motif_env
# Load shared functions
source ~/GitHub/Barley_Mutated/02_analysis/mutation_motif/prep_counts_functions.sh

# User provided input arguments
VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/hybrid13.SNPs.common.vcf.gz"
# Reference fasta file
REF_FASTA="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta"
REF_FAI="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta.fai"
# Size around SNV to expand window in bp (min, 1)
#   Note: Integer will be added in each direction
#   Tip: If you run aln_to_counts and get the error: "ValueError: not all sequences have same length"
#       keep reducing the WIN_FLANK_SIZE until that error goes away
WIN_FLANK_SIZE="100"
# The number of bases per side to include
FLANK_SIZE="2"
# Where should we output our files?
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/mutation_motif/hybrid13.SNPs.common"
# Temporary directory for intermediate files
TEMP_DIR="/scratch.global/liux1299/temp_mut_motif.hybrid13.SNPs.common"

#----------------------
# Create output directories and subdirectories
mkdir -p ${OUT_DIR} ${OUT_DIR}/counts_tables
mkdir -p ${TEMP_DIR} ${TEMP_DIR}/split_vcf ${TEMP_DIR}/split_bed ${TEMP_DIR}/split_fasta

# Prepare output prefix of VCF file
if [[ "${VCF}" == *"gz"* ]]; then
    # We are working with gzipped VCF file
    out_prefix=$(basename ${VCF} .vcf.gz)
else
    # We are working with uncompressed vcf
    out_prefix=$(basename ${VCF} .vcf)
fi

# Split VCF file based on mutation motif categories for the aln_to_counts option:
#   --direction [AtoC|AtoG|AtoT|CtoA|CtoG|CtoT|GtoA|GtoC|GtoT|TtoA|TtoC|TtoG]
bcftools view -i 'REF="A" & ALT="C"' ${VCF} -O v -o ${TEMP_DIR}/split_vcf/AtoC_${out_prefix}.vcf
bcftools view -i 'REF="A" & ALT="G"' ${VCF} -O v -o ${TEMP_DIR}/split_vcf/AtoG_${out_prefix}.vcf
bcftools view -i 'REF="A" & ALT="T"' ${VCF} -O v -o ${TEMP_DIR}/split_vcf/AtoT_${out_prefix}.vcf
bcftools view -i 'REF="C" & ALT="A"' ${VCF} -O v -o ${TEMP_DIR}/split_vcf/CtoA_${out_prefix}.vcf
bcftools view -i 'REF="C" & ALT="G"' ${VCF} -O v -o ${TEMP_DIR}/split_vcf/CtoG_${out_prefix}.vcf
bcftools view -i 'REF="C" & ALT="T"' ${VCF} -O v -o ${TEMP_DIR}/split_vcf/CtoT_${out_prefix}.vcf
bcftools view -i 'REF="G" & ALT="A"' ${VCF} -O v -o ${TEMP_DIR}/split_vcf/GtoA_${out_prefix}.vcf
bcftools view -i 'REF="G" & ALT="C"' ${VCF} -O v -o ${TEMP_DIR}/split_vcf/GtoC_${out_prefix}.vcf
bcftools view -i 'REF="G" & ALT="T"' ${VCF} -O v -o ${TEMP_DIR}/split_vcf/GtoT_${out_prefix}.vcf
bcftools view -i 'REF="T" & ALT="A"' ${VCF} -O v -o ${TEMP_DIR}/split_vcf/TtoA_${out_prefix}.vcf
bcftools view -i 'REF="T" & ALT="C"' ${VCF} -O v -o ${TEMP_DIR}/split_vcf/TtoC_${out_prefix}.vcf
bcftools view -i 'REF="T" & ALT="G"' ${VCF} -O v -o ${TEMP_DIR}/split_vcf/TtoG_${out_prefix}.vcf

# Prepare array of split VCF files
SPLIT_VCF_ARR=($(realpath ${TEMP_DIR}/split_vcf/*${out_prefix}*))

# Convert split vcf files to BED format and extend interval by WIN_FLANK_SIZE in each direction
parallel --verbose vcf_to_bed {} ${TEMP_DIR}/split_bed ${REF_FAI} ${WIN_FLANK_SIZE} ::: ${SPLIT_VCF_ARR[@]}
# Prepare array of extended interval BED files
SPLIT_BED_ARR=($(realpath ${TEMP_DIR}/split_bed/*.winFlank${WIN_FLANK_SIZE}bp.bed))

# Extract sequences from fasta as defined in BED
parallel --verbose bed_to_fasta {} ${REF_FASTA} ${TEMP_DIR}/split_fasta ::: ${SPLIT_BED_ARR[@]}
# Prepare array of FASTA files
SPLIT_FASTA_ARR=($(realpath ${TEMP_DIR}/split_fasta/*.winFlank${WIN_FLANK_SIZE}bp.fasta))

# Prepare counts table format
parallel --verbose run_aln_to_counts {} ${OUT_DIR}/counts_tables ${FLANK_SIZE} ${out_prefix} ${WIN_FLANK_SIZE} ::: ${SPLIT_FASTA_ARR[@]}

# Combine separate counts table into single file that is suitable for spectra analysis
all_counts \
    --counts_pattern "${OUT_DIR}/counts_tables/*${out_prefix}.winFlank${WIN_FLANK_SIZE}bp.txt*" \
    --output_path ${OUT_DIR}/counts_tables \
    --strand_symmetric \
    --force_overwrite
# Rename combined counts output file
mv ${OUT_DIR}/counts_tables/combined_counts.txt ${OUT_DIR}/counts_tables/combined_counts_${out_prefix}.winFlank${WIN_FLANK_SIZE}bp.txt
mv ${OUT_DIR}/counts_tables/combined_counts.log ${OUT_DIR}/counts_tables/combined_counts_${out_prefix}.winFlank${WIN_FLANK_SIZE}bp.log
