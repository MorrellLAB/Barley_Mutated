# Compare VeP premature stop to Annovar

library(purrr)
library(dplyr)
library(tidyverse)

# User provided input arguments
vep_snps_fp <- "~/Dropbox/Projects/Barley_Mutated/analyses/VEP/HC_LC_gff_SNPs_private_all_samples/mut8_and_3mut10xGenomics.SNPs.private.txt"
vep_indels_fp <- "~/Dropbox/Projects/Barley_Mutated/analyses/VEP/HC_LC_gff_INDELs_private_all_samples/mut8_and_3mut10xGenomics.INDELs.private.txt"

annovar_snps_fp <- "~/Dropbox/Projects/Barley_Mutated/analyses/Annovar/all/mut8_and_3mut10xGenomics.SNPs.private_annovar_input.txt.exonic_variant_function"
annovar_indels_fp <- "~/Dropbox/Projects/Barley_Mutated/analyses/Annovar/all/mut8_and_3mut10xGenomics.INDELs.private_annovar_input.txt.exonic_variant_function"

#---------------
# Read in data
# Only keep unique rows
vep_snps_df <- read.delim(file=vep_snps_fp, header=FALSE, sep="\t", comment.char="#") %>% distinct()
# Add column names
# Formatted column names in BASH using:
# grep "#Uploaded_variation" mut8_and_3mut10xGenomics.SNPs.private.txt | tr '\t' '\n' | sed 's/#//' | sed 's/^/"/' | sed 's/$/"/' | tr '\n' ',' | sed 's/,/, /g'
colnames(vep_snps_df) <- c("Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "Extra")

vep_indels_df <- read.delim(file=vep_indels_fp, header=FALSE, sep="\t", comment.char="#") %>% distinct()
# Add column names
colnames(vep_indels_df) <- c("Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "Extra")

# Reformat VeP data for easier comparison
# SNPs
vep_split_chr_pos <- data.frame(do.call("rbind", strsplit(as.character(vep_snps_df$Location), ":", fix=TRUE)))
colnames(vep_split_chr_pos) <- c("chr", "pos")
vep_snps_df <- cbind(vep_split_chr_pos, vep_snps_df)
# indels
vep_split_chr_pos_indels <- data.frame(do.call("rbind", strsplit(as.character(vep_indels_df$Location), ":", fix=TRUE)))
colnames(vep_split_chr_pos_indels) <- c("chr", "pos")
vep_pos_split_indels <- data.frame(do.call("rbind", strsplit(as.character(vep_split_chr_pos_indels$pos), "-", fix=TRUE)))
colnames(vep_pos_split_indels) <- c("start_pos", "end_pos")
vep_indels_df <- cbind(chr=vep_split_chr_pos_indels$chr, vep_pos_split_indels, vep_indels_df)

# Read in Annovar files
annovar_snps_df <- read.delim(file=annovar_snps_fp, header=FALSE, sep="\t", comment.char="#")
annovar_indels_df <- read.delim(file=annovar_indels_fp, header=FALSE, sep="\t", comment.char="#")
# Reformat annovar columns
# SNPs
annovar_snps_split_feature <- data.frame(as.character(map(strsplit(annovar_snps_df$V3, split = ":"), 2)))
colnames(annovar_snps_split_feature) <- "Annovar_Feature"
annovar_snps_df <- cbind(annovar_snps_split_feature, annovar_snps_df)
# indels
annovar_indels_split_feature <- data.frame(as.character(map(strsplit(annovar_indels_df$V3, split = ":"), 2)))
colnames(annovar_indels_split_feature) <- "Annovar_Feature"
annovar_indels_df <- cbind(annovar_indels_split_feature, annovar_indels_df)

# Add data source labels
vep_snps_df$data_source <- "VeP"
vep_indels_df$data_source <- "VeP"
annovar_snps_df$data_source <- "Annovar"
annovar_indels_df$data_source <- "Annovar"

# Merge data frames and bring data source columns to front
# Remove intergenic variants for comparison purposes
snp_table <- merge(vep_snps_df, annovar_snps_df, by.x=c("Feature", "chr", "pos"), by.y=c("Annovar_Feature", "V4", "V5"), all.x=TRUE, all.y=TRUE) %>%
  select(data_source.x, data_source.y, everything()) %>%
  filter(Consequence!="intergenic_variant")

indel_table <- merge(vep_indels_df, annovar_indels_df, by.x=c("chr", "start_pos", "end_pos"), by.y=c("V4", "V5", "V6"), all.x=TRUE, all.y=TRUE) %>%
  select(data_source.x, data_source.y, everything()) %>%
  filter(Consequence!="intergenic_variant")

# SNPs - look at stop gains
# Move relevant columns to front
stopgain_snp_table <- snp_table %>%
  filter(Consequence=="stop_gained" | V2=="stopgain") %>%
  select(data_source.x, Consequence, data_source.y, V2, everything())
# Check if Annovar calls stopgain but VeP doesn't
snp_table %>%
  filter(Consequence!="stop_gained" & V2=="stopgain") %>%
  select(data_source.x, Consequence, data_source.y, V2, everything())

# Check indels
# Look at positions that Annovar called (since we're only looking at exonic_variants)
annovar_called_indels <- indel_table %>%
  filter(data_source.y=="Annovar") %>%
  select(data_source.x, Consequence, data_source.y, V2, everything())
