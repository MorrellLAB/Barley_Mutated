#!/usr/bin/env Rscript

# Generates summary plot of VEP result counts.
# This script is for initial exploration purposes.

# Dependencies
library(tidyverse)
library(scales)
library(RColorBrewer)
library(ggthemes)
require(gridExtra)
require(grid)

# User defined input arguments
m01_fp <- "~/Dropbox/Projects/Barley_Mutated/analyses/VEP/per_sample-singletons/M01-3-3_filtered_singletons_only.txt"
m20_fp <- "~/Dropbox/Projects/Barley_Mutated/analyses/VEP/per_sample-singletons/M20-2-2_filtered_singletons_only.txt"
m29_fp <- "~/Dropbox/Projects/Barley_Mutated/analyses/VEP/per_sample-singletons/M29-2-2_filtered_singletons_only.txt"
# Where do we want to save our plots?
out_dir <- "~/Dropbox/Projects/Barley_Mutated/analyses/VEP/plots"

# Generate output plot filenames
tc.outfilename <- paste0(out_dir, "/vep_total_count_singletons.pdf")
psc.outfilename <- paste0(out_dir, "/vep_per_sample_count_singletons.pdf")

#---------------------------
# Prepare headers
# First 30 lines start with "##" + 1 line starts with #Uploaded_variation
num_header_lines <- 31
custom_header <- c("Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "Extra")

# Read in data
m01_df <- read.csv(m01_fp, sep="\t", skip=num_header_lines, header=FALSE, col.names = custom_header)
m20_df <- read.csv(m20_fp, sep="\t", skip=num_header_lines, header=FALSE, col.names = custom_header)
m29_df <- read.csv(m29_fp, sep="\t", skip=num_header_lines, header=FALSE, col.names = custom_header)

# Add column to distinguish each sample
m01_df["sample_id"] <- "M01-3-3"
m20_df["sample_id"] <- "M20-2-2"
m29_df["sample_id"] <- "M29-2-2"

# Combine data into single data frame for easier plotting
df <- rbind(m01_df, m20_df, m29_df)
# Subset and remove largest class (intergenic variants)
df_noIntrg <- df[df$Consequence != "intergenic_variant", ]
# Subset and remove largest 3 (top 3) classes:
#   1) intergenic variants
#   2) downstream_gene_variant
#   3) upstream_gene_variant
df_minusTop3 <- filter(df, Consequence != "intergenic_variant" & Consequence != "downstream_gene_variant" & Consequence != "upstream_gene_variant")

# Count number of each consequence type
# Total counts of all samples combined
total_counts_df <- df %>% 
  group_by(Consequence) %>% 
  summarise(counts = n()) %>%
  mutate(percentage = signif(prop.table(counts), digits=2))

# Total counts, excluding intergenic variants
total_counts_noIntrg <- df_noIntrg %>%
  group_by(Consequence) %>% 
  summarise(counts = n()) %>%
  mutate(percentage = signif(prop.table(counts), digits=2))

# Total counts, excluding top 3 largest classes
total_counts_minusTop3 <- df_minusTop3 %>%
  group_by(Consequence) %>% 
  summarise(counts = n()) %>%
  mutate(percentage = signif(prop.table(counts), digits=2))

# Counts per sample
# Note: proportion is relative to each sample (i.e., consequence_type_sample/sample_total)
samp_counts_df <- df %>% 
  group_by(Consequence, sample_id) %>% 
  summarise(counts_per_sample = n()) %>%
  group_by(sample_id) %>%
  mutate(percentage = signif(prop.table(counts_per_sample), digits=2))

# Counts per sample, excluding intergenic variants
samp_counts_noIntrg <- df_noIntrg %>%
  group_by(Consequence, sample_id) %>% 
  summarise(counts_per_sample = n()) %>%
  group_by(sample_id) %>%
  mutate(percentage = signif(prop.table(counts_per_sample), digits=2))

# Counts per sample, excluding top 3 largest classes
samp_counts_minusTop3 <- df_minusTop3 %>%
  group_by(Consequence, sample_id) %>% 
  summarise(counts_per_sample = n()) %>%
  group_by(sample_id) %>%
  mutate(percentage = signif(prop.table(counts_per_sample), digits=2))

################################
# Generate the plots
# Shared plot options
plot_options <- list(theme_minimal(),
                     ylab("Consequence Type"),
                     xlab("Count"),
                     theme(plot.title=element_text(hjust=0.5, size=14, face="bold"),
                           axis.text=element_text(size=12),
                           axis.title=element_text(size=14),
                           # Legend
                           legend.text=element_text(size=11))
                     )

# Total counts, all samples combined and all consequence types included
# Counts label
tc.all <- ggplot(data=total_counts_df, aes(x=counts, y=reorder(Consequence, counts), label=counts)) +
  geom_col() +
  geom_text(hjust=-0.2) +
  scale_x_continuous(breaks=seq(0, 300000, by=50000)) +
  ggtitle("All samples combined with all consequence types") +
  plot_options
# Preview plot
tc.all

# Percentage label
# ggplot(data=total_counts_df, aes(x=counts, y=reorder(Consequence, -counts), label=scales::percent(percentage))) +
#   geom_col() +
#   geom_text(hjust=-0.2) +
#   xlim(0, 270000)

# Total counts, exlcude intergenic variants
tc.noIntrg <- ggplot(data=total_counts_noIntrg, aes(x=counts, y=reorder(Consequence, counts), label=counts)) +
  geom_col() +
  geom_text(hjust=-0.2) +
  plot_options
# Preview plot
tc.noIntrg

# Total counts, exlcude top 3 largest classes
tc.minusTop3 <- ggplot(data=total_counts_minusTop3, aes(x=counts, y=reorder(Consequence, counts), label=counts)) +
  geom_col() +
  geom_text(hjust=-0.2) +
  ggtitle("All samples combined excluding 3 largest consequence classes") +
  plot_options
# Preview plot
tc.minusTop3

# Per sample counts, all consequence types included
# Counts label
psc.all <- ggplot(data=samp_counts_df, aes(x=counts_per_sample, y=reorder(Consequence, counts_per_sample), fill=sample_id, group=sample_id, label=counts_per_sample)) +
  geom_col(position="dodge") +
  geom_text(position=position_dodge(width=0.9), hjust=-0.2) +
  scale_fill_viridis_d(option = "D") +
  guides(fill=guide_legend(title="SampleID")) +
  ggtitle("Counts per sample with all consequence types") +
  plot_options
# Preview plot
psc.all

# Per sample counts, exclude intergenic variants
psc.noIntrg <- ggplot(data=samp_counts_noIntrg, aes(x=counts_per_sample, y=reorder(Consequence, counts_per_sample), fill=sample_id, group=sample_id, label=counts_per_sample)) +
  geom_col(position="dodge") +
  geom_text(position=position_dodge(width=0.9), hjust=-0.2) +
  scale_fill_viridis_d(option = "D") +
  guides(fill=guide_legend(title="SampleID")) +
  plot_options
# Preview plot
psc.noIntrg

# Per sample counts, exclude intergenic variants
psc.minusTop3 <- ggplot(data=samp_counts_minusTop3, aes(x=counts_per_sample, y=reorder(Consequence, counts_per_sample), fill=sample_id, group=sample_id, label=counts_per_sample)) +
  geom_col(position="dodge") +
  geom_text(position=position_dodge(width=0.9), hjust=-0.2) +
  scale_fill_viridis_d(option = "D") +
  guides(fill=guide_legend(title="SampleID")) +
  ggtitle("Counts per sample excluding 3 largest consequence classes") +
  plot_options
# Preview plot
psc.minusTop3

####################
# Save plot drafts to file
tc.outplot <- arrangeGrob(tc.all, tc.minusTop3, nrow=2)
ggsave(filename=tc.outfilename, tc.outplot, width = 20, height = 12)

psc.outplot <- arrangeGrob(psc.all, psc.minusTop3, nrow=2)
ggsave(filename=psc.outfilename, psc.outplot, width = 20, height = 14)
