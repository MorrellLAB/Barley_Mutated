#!/usr/bin/env Rscript

library(tidyverse)
library(tidyr)
library(ggplot2)
library(ggpubr)

# Set directory to output plots
setwd("~/Dropbox/Projects/Barley_Mutated/analyses/yield_testing")

# Filepaths
# Plot level averages
fp_plot_level_avg = "~/Dropbox/Projects/Barley_Mutated/field_data_and_planting/Field_Testing/Raw_and_Spatially_Adjusted_2020-2022/mmx 2020-2022_plot level_031623.csv"
# Field level averages
fp_trial_avg = "~/Dropbox/Projects/Barley_Mutated/field_data_and_planting/Field_Testing/Raw_and_Spatially_Adjusted_2020-2022/mmx 2020-2022_trial avg_031623.csv"
dsnps_counts_fp = "~/Dropbox/Projects/Barley_Mutated/analyses/bad_mutations/visualization/mut_private_per_sample_counts.del_vs_tol.nonsyn_vs_syn.txt"

out_dir = "~/Dropbox/Projects/Barley_Mutated/analyses/bad_mutations/plots"
out_prefix = "corr_mut_private_per_sample_counts"

# Check varieties and Morex parent
check_lines <- c("CONLON", "FEG141-20", "LACEY", "ND_GENESIS", "ND20448", "ND26104", "PINNACLE", "RASMUSSON")

# Sequenced mutated lines
seq_mut_lines <- c("M01-3-3", "M02-1", "M11-3-3", "M14-2-2", "M20-2-2", "M28-3-1", "M29-2-2", "M35-3-2", "M36-1-2", "M39-2-3", "M41-2-1")

#-----------------
# Trial level average
df <- read.csv(file=fp_trial_avg, header=TRUE, skip=1)
# Remove "FILLER" rows
df <- df %>% filter(!grepl("FILLER", line_name))

# Convert yield to numeric
df$avg_yield_bua_use <- as.numeric(df$avg_yield_bua_use)
df$avg_yield_kgha_use <- as.numeric(df$avg_yield_kgha_use)
df$avg_yield_kgha <- as.numeric(df$avg_yield_kgha)

# Highlight mutated lines
df <- df %>% mutate(linetype=case_when((line_name == "MOREX_W2017") ~ "Parent",
                                       (line_name %in% check_lines) ~ "Check",
                                       !(line_name %in% check_lines) ~ "Sodium azide"))

# Shade sequenced mutated lines
df <- df %>% mutate(seq_status=ifelse(line_name %in% seq_mut_lines, "Sequenced", "Not Sequenced"))

# Plot level average
# Read dsnps counts table
dsnps_counts <- read.delim(dsnps_counts_fp, header = TRUE)
# Create column of short names
dsnps_counts_short_names <- separate(dsnps_counts, col=sample, into=c("short_name"), sep="-")
# Add short names to original table
dsnps_counts$short_name <- dsnps_counts_short_names$short_name
dsnps_counts <- dsnps_counts %>% select(short_name, everything())

# Plot level average
df.pl <- read.csv(file=fp_plot_level_avg, header=TRUE, skip=1)
# Remove "FILLER" rows
df.pl <- df.pl %>% filter(!grepl("FILLER", line_name))

# Convert yield to numeric
df.pl$yield_bua_use <- as.numeric(df.pl$yield_bua_use)
df.pl$yield_kgha_use <- as.numeric(df.pl$yield_kgha_use)
df.pl$yield_kgha <- as.numeric(df.pl$yield_kgha)
# Convert heading_dap to numeric
df.pl$heading_dap_use <- as.numeric(df.pl$heading_dap_use)
# Convert height to numeric
df.pl$height_cm_use <- as.numeric(df.pl$height_cm_use)

# Highlight mutated lines
df.pl <- df.pl %>% mutate(linetype=case_when((line_name == "MOREX_W2017") ~ "Parent",
                                             (line_name %in% check_lines) ~ "Check",
                                             !(line_name %in% check_lines) ~ "Sodium azide"))

# Shade sequenced mutated lines
df.pl <- df.pl %>% mutate(seq_status=ifelse(line_name %in% seq_mut_lines, "Sequenced", "Not Sequenced"))

# Create column of short names
uniq_line_names <- sort(unique(df.pl$line_name))
short_names_df <- data.frame(short_name=sort(vapply(strsplit(uniq_line_names, split="-"), '[', 1, FUN.VALUE=character(1))), line_name=uniq_line_names)

############
# Trait means
tmp.df.pl.yield <- df.pl %>%
  select(c(line_name, yield_kgha_use, yield_bua_use, yield_kgha)) %>%
  group_by(line_name) %>%
  select(where(is.numeric)) %>%
  summarise(across(everything(), mean))

# Remove rows that contain NA
tmp.df.pl.hdap <- df.pl %>%
  select(c(line_name, heading_dap_use)) %>%
  filter(complete.cases(.)) %>%
  group_by(line_name) %>%
  summarise(mean_heading_dap_use=(mean(heading_dap_use)))

tmp.df.pl.height <- df.pl %>%
  select(c(line_name, height_cm_use)) %>%
  filter(complete.cases(.)) %>%
  group_by(line_name) %>%
  summarise(mean_height_cm_use=(mean(height_cm_use)))

# Combine trait means into single df
# And add short names
df.pl.trait_means <- left_join(short_names_df, left_join(tmp.df.pl.yield, left_join(tmp.df.pl.hdap, tmp.df.pl.height, by="line_name"), by="line_name"), by="line_name")

############
# Plot level averages
# Per sample mean and move to first column
dsnps_trait_means <- left_join(x=dsnps_counts, y=df.pl.trait_means, by="short_name") %>%
  select(short_name, everything())

# Test for normality - all traits
# Traits
shapiro.test(dsnps_trait_means$yield_kgha_use)
shapiro.test(dsnps_trait_means$yield_bua_use)
shapiro.test(dsnps_trait_means$yield_kgha)
shapiro.test(dsnps_trait_means$mean_heading_dap_use)
shapiro.test(dsnps_trait_means$mean_height_cm_use)
# Counts
shapiro.test(dsnps_trait_means$num_snps)
shapiro.test(dsnps_trait_means$num_hom_alt)
# p-value = 0.01101
shapiro.test(dsnps_trait_means$num_het)
# p-value = 0.04536
shapiro.test(dsnps_trait_means$total_del)
shapiro.test(dsnps_trait_means$total_tol)
shapiro.test(dsnps_trait_means$total_syn)
shapiro.test(dsnps_trait_means$total_nonsyn)
# Check QQ plots - all traits
# Traits
ggqqplot(dsnps_trait_means$yield_kgha_use, ylab="Mean Grain Yield (kg/ha)")
ggqqplot(dsnps_trait_means$mean_heading_dap_use, ylab="Mean Heading DAP")
ggqqplot(dsnps_trait_means$mean_height_cm_use, ylab="Mean Height (cm)")
# Counts
ggqqplot(dsnps_trait_means$num_snps, ylab="# SNPs")
ggqqplot(dsnps_trait_means$num_het, ylab="# HET SNPs")
ggqqplot(dsnps_trait_means$total_del, ylab="# deleterious SNPs")
ggqqplot(dsnps_trait_means$total_tol, ylab="# tolerated SNPs")
ggqqplot(dsnps_trait_means$total_syn, ylab="# synonynmous SNPs")
ggqqplot(dsnps_trait_means$total_nonsyn, ylab="# nonsynonynmous SNPs")

# Correlation test
res_del_vs_yield <- cor.test(dsnps_trait_means$total_del, dsnps_trait_means$yield_kgha_use, method="spearman", alternative="two.sided")
res_del_vs_yield
res_del_vs_yield$p.value

# Prepare base filepath
out_fp_base = paste0(out_dir, '/', out_prefix)

# Number of deleterious
jpeg(filename=paste0(out_fp_base, '-total_del_vs_yield_kgha_use.jpg'), width=8, height=6, units='in', res=300)
ggscatter(dsnps_trait_means, x="total_del", y="yield_kgha_use",
          add="reg.line", conf.int=TRUE,
          cor.coef=TRUE, cor.method="spearman",
          xlab="Number of deleterious SNPs", ylab="Yield (Kg/ha)")
dev.off()

jpeg(filename=paste0(out_fp_base, '-total_del_vs_mean_heading_dap_use.jpg'), width=8, height=6, units='in', res=300)
ggscatter(dsnps_trait_means, x="total_del", y="mean_heading_dap_use",
          add="reg.line", conf.int=TRUE,
          cor.coef=TRUE, cor.method="spearman",
          xlab="Number of deleterious SNPs", ylab="Mean Heading Days After Planting")
dev.off()

jpeg(filename=paste0(out_fp_base, '-total_del_vs_mean_height_cm_use.jpg'), width=8, height=6, units='in', res=300)
ggscatter(dsnps_trait_means, x="total_del", y="mean_height_cm_use",
          add="reg.line", conf.int=TRUE,
          cor.coef=TRUE, cor.method="spearman",
          xlab="Number of deleterious SNPs", ylab="Mean Height (cm)")
dev.off()

# Number of nonsynonymous
jpeg(filename=paste0(out_fp_base, '-total_nonsyn_vs_yield_kgha_use.jpg'), width=8, height=6, units='in', res=300)
ggscatter(dsnps_trait_means, x="total_nonsyn", y="yield_kgha_use",
          add="reg.line", conf.int=TRUE,
          cor.coef=TRUE, cor.method="pearson",
          xlab="Number of nonsynonymous SNPs", ylab="Yield (Kg/ha)")
dev.off()

jpeg(filename=paste0(out_fp_base, '-total_nonsyn_vs_mean_heading_dap_use.jpg'), width=8, height=6, units='in', res=300)
ggscatter(dsnps_trait_means, x="total_nonsyn", y="mean_heading_dap_use",
          add="reg.line", conf.int=TRUE,
          cor.coef=TRUE, cor.method="pearson",
          xlab="Number of nonsynonymous SNPs", ylab="Mean Heading Days After Planting")
dev.off()

jpeg(filename=paste0(out_fp_base, '-total_nonsyn_vs_mean_height_cm_use.jpg'), width=8, height=6, units='in', res=300)
ggscatter(dsnps_trait_means, x="total_nonsyn", y="mean_height_cm_use",
          add="reg.line", conf.int=TRUE,
          cor.coef=TRUE, cor.method="pearson",
          xlab="Number of nonsynonymous SNPs", ylab="Mean Height (cm)")
dev.off()

# Number of tolerated
jpeg(filename=paste0(out_fp_base, '-total_tol_vs_yield_kgha_use.jpg'), width=8, height=6, units='in', res=300)
ggscatter(dsnps_trait_means, x="total_tol", y="yield_kgha_use",
          add="reg.line", conf.int=TRUE,
          cor.coef=TRUE, cor.method="pearson",
          xlab="Number of tolerated SNPs", ylab="Yield (Kg/ha)")
dev.off()

jpeg(filename=paste0(out_fp_base, '-total_tol_vs_mean_heading_dap_use.jpg'), width=8, height=6, units='in', res=300)
ggscatter(dsnps_trait_means, x="total_tol", y="mean_heading_dap_use",
          add="reg.line", conf.int=TRUE,
          cor.coef=TRUE, cor.method="pearson",
          xlab="Number of tolerated SNPs", ylab="Mean Heading Days After Planting")
dev.off()

jpeg(filename=paste0(out_fp_base, '-total_tol_vs_mean_height_cm_use.jpg'), width=8, height=6, units='in', res=300)
ggscatter(dsnps_trait_means, x="total_tol", y="mean_height_cm_use",
          add="reg.line", conf.int=TRUE,
          cor.coef=TRUE, cor.method="pearson",
          xlab="Number of tolerated SNPs", ylab="Mean Height (cm)")
dev.off()

# Number of SNPs
jpeg(filename=paste0(out_fp_base, '-snps_vs_yield_kgha_use.jpg'), width=8, height=6, units='in', res=300)
ggscatter(dsnps_trait_means, x="num_snps", y="yield_kgha_use",
          add="reg.line", conf.int=TRUE,
          cor.coef=TRUE, cor.method="pearson",
          xlab="Number of SNPs", ylab="Yield (Kg/ha)")
dev.off()

jpeg(filename=paste0(out_fp_base, '-snps_vs_mean_heading_dap_use.jpg'), width=8, height=6, units='in', res=300)
ggscatter(dsnps_trait_means, x="num_snps", y="mean_heading_dap_use",
          add="reg.line", conf.int=TRUE,
          cor.coef=TRUE, cor.method="pearson",
          xlab="Number of SNPs", ylab="Mean Heading Days After Planting")
dev.off()

jpeg(filename=paste0(out_fp_base, '-snps_vs_mean_height_cm_use.jpg'), width=8, height=6, units='in', res=300)
ggscatter(dsnps_trait_means, x="num_snps", y="mean_height_cm_use",
          add="reg.line", conf.int=TRUE,
          cor.coef=TRUE, cor.method="pearson",
          xlab="Number of SNPs", ylab="Mean Height (cm)")
dev.off()

# Number of synonymous
jpeg(filename=paste0(out_fp_base, '-total_syn_vs_yield_kgha_use.jpg'), width=8, height=6, units='in', res=300)
ggscatter(dsnps_trait_means, x="total_syn", y="yield_kgha_use",
          add="reg.line", conf.int=TRUE,
          cor.coef=TRUE, cor.method="pearson",
          xlab="Number of synonymous SNPs", ylab="Yield (Kg/ha)")
dev.off()

jpeg(filename=paste0(out_fp_base, '-total_syn_vs_mean_heading_dap_use.jpg'), width=8, height=6, units='in', res=300)
ggscatter(dsnps_trait_means, x="total_syn", y="mean_heading_dap_use",
          add="reg.line", conf.int=TRUE,
          cor.coef=TRUE, cor.method="pearson",
          xlab="Number of synonymous SNPs", ylab="Mean Heading Days After Planting")
dev.off()

jpeg(filename=paste0(out_fp_base, '-total_syn_vs_mean_height_cm_use.jpg'), width=8, height=6, units='in', res=300)
ggscatter(dsnps_trait_means, x="total_syn", y="mean_height_cm_use",
          add="reg.line", conf.int=TRUE,
          cor.coef=TRUE, cor.method="pearson",
          xlab="Number of synonymous SNPs", ylab="Mean Height (cm)")
dev.off()
