#!/usr/bin/env Rscript

library(dplyr)

args = commandArgs(trailingOnly=TRUE)

# Merge BAD_Mutations compiled predict report and VeP report
# We'll primarily need the amino acid info from VeP, but we'll keep the other info
# for other downstream uses

# User provided input arguments
bad_mut_predict_fp <- args[1]
vep_fp <- args[2]
out_fp <- args[3]

#----------------------
# Read in data
predict_df <- read.delim(file=bad_mut_predict_fp, header=TRUE, sep="\t")
# Load VeP report and keep unique rows
vep_df <- read.table(file=vep_fp, header=FALSE, sep="\t", stringsAsFactors=FALSE, na.strings=c("NA", '-')) %>% distinct()
colnames(vep_df) <- c("VariantID", "Position", "Allele", "GeneID", "Feature", "Feature_type", "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "Extra")

# Only include rows that exist in both BAD_Mutations report and VeP report
snp_table <- merge(predict_df, vep_df, by.x=c("VariantID", "GeneID"), by.y=c("VariantID", "Feature"))

# Print to stdout
write.table(x=snp_table, file=out_fp, quote=FALSE, sep="\t", row.names=FALSE)
