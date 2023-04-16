#!/usr/bin/env Rscript

library(tidyverse)
library(grid)
library(gridExtra)
library(lattice)

# Set directory to output plots
setwd("~/Dropbox/Projects/Barley_Mutated/analyses/yield_testing")

# Filepaths
# Field level averages
fp_trial_avg = "~/Dropbox/Projects/Barley_Mutated/field_data_and_planting/Field_Testing/Raw_and_Spatially_Adjusted_2020-2022/mmx 2020-2022_trial avg_031623.csv"
# Plot level averages
fp_plot_level_avg = "~/Dropbox/Projects/Barley_Mutated/field_data_and_planting/Field_Testing/Raw_and_Spatially_Adjusted_2020-2022/mmx 2020-2022_plot level_031623.csv"

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

### All locations and years combined
# Spatially adjusted
p1.sa.avg_yield <- ggplot(data=df, aes(x=reorder(line_name, avg_yield_kgha_use, FUN=median), y=as.numeric(avg_yield_kgha_use), color=linetype, fill=seq_status, alpha=seq_status)) +
  geom_boxplot(width=0.5, outlier.shape=NA) +
  geom_jitter(width=0, size=1, alpha=0.8) +
  scale_color_manual(values=c("grey80", "midnightblue", "darkred")) +
  scale_fill_manual(values=c("grey80", "darkred")) +
  scale_alpha_manual(values=c(0.001, 0.4)) +
  xlab("Line name") +
  ylab("Average Grain Yield (Kg/ha)") +
  ggtitle("Trial Average - Spatially Adjusted") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        plot.title=element_text(hjust=0.5))

# Raw data (without spatial adjustment)
p2.raw.avg_yield <- ggplot(data=df, aes(x=reorder(line_name, avg_yield_kgha, FUN=median), y=as.numeric(avg_yield_kgha), color=linetype, fill=seq_status, alpha=seq_status)) +
  geom_boxplot(width=0.5, outlier.shape=NA) +
  geom_jitter(width=0, size=1, alpha=0.8) +
  scale_color_manual(values=c("grey80", "midnightblue", "darkred")) +
  scale_fill_manual(values=c("grey80", "darkred")) +
  scale_alpha_manual(values=c(0.001, 0.4)) +
  xlab("Line name") +
  ylab("Average Grain Yield (Kg/ha)") +
  ggtitle("Trial Average - Raw") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        plot.title=element_text(hjust=0.5))

grid.arrange(p1.sa.avg_yield, p2.raw.avg_yield, nrow=2)
combo.avg_yield <- arrangeGrob(p1.sa.avg_yield, p2.raw.avg_yield, nrow=2)
# Save plot
ggsave("trial_avg_spatially_adjusted_vs_raw.png", combo.avg_yield, dpi=300)

#-----------------
# Plot level average
df.pl <- read.csv(file=fp_plot_level_avg, header=TRUE, skip=1)
# Remove "FILLER" rows
df.pl <- df.pl %>% filter(!grepl("FILLER", line_name))

# Convert yield to numeric
df.pl$yield_bua_use <- as.numeric(df.pl$yield_bua_use)
df.pl$yield_kgha_use <- as.numeric(df.pl$yield_kgha_use)
df.pl$yield_kgha <- as.numeric(df.pl$yield_kgha)

# Highlight mutated lines
df.pl <- df.pl %>% mutate(linetype=case_when((line_name == "MOREX_W2017") ~ "Parent",
                                         (line_name %in% check_lines) ~ "Check",
                                         !(line_name %in% check_lines) ~ "Sodium azide"))

# Shade sequenced mutated lines
df.pl <- df.pl %>% mutate(seq_status=ifelse(line_name %in% seq_mut_lines, "Sequenced", "Not Sequenced"))

### All locations and years combined
# Spatially adjusted
p1.sa.yield <- ggplot(data=df.pl, aes(x=reorder(line_name, yield_kgha_use, FUN=median), y=as.numeric(yield_kgha_use), color=linetype, fill=seq_status, alpha=seq_status)) +
  geom_boxplot(width=0.5, outlier.shape=NA) +
  geom_jitter(width=0.2, size=1, alpha=0.8) +
  scale_color_manual(values=c("grey80", "midnightblue", "darkred")) +
  scale_fill_manual(values=c("grey80", "darkred")) +
  scale_alpha_manual(values=c(0.001, 0.4)) +
  xlab("Line name") +
  ylab("Grain Yield (Kg/ha)") +
  ggtitle("Plot Level - Spatially Adjusted") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        plot.title=element_text(hjust=0.5))

# Raw data (without spatial adjustment)
p2.raw.yield <- ggplot(data=df.pl, aes(x=reorder(line_name, yield_kgha, FUN=median), y=as.numeric(yield_kgha), color=linetype, fill=seq_status, alpha=seq_status)) +
  geom_boxplot(width=0.5, outlier.shape=NA) +
  geom_jitter(width=0.2, size=1, alpha=0.8) +
  scale_color_manual(values=c("grey80", "midnightblue", "darkred")) +
  scale_fill_manual(values=c("grey80", "darkred")) +
  scale_alpha_manual(values=c(0.001, 0.4)) +
  xlab("Line name") +
  ylab("Grain Yield (Kg/ha)") +
  ggtitle("Plot Level - Raw") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        plot.title=element_text(hjust=0.5))

grid.arrange(p1.sa.yield, p2.raw.yield, nrow=2)
combo.plot_yield <- arrangeGrob(p1.sa.yield, p2.raw.yield, nrow=2)
# Save plot
ggsave("plot_level_spatially_adjusted_vs_raw.png", combo.plot_yield, dpi=300)

### Plot level by year
ggplot(data=df.pl, aes(x=reorder(line_name, yield_kgha_use, FUN=median), y=as.numeric(yield_kgha_use), color=linetype, fill=seq_status, alpha=seq_status)) +
  geom_boxplot(width=0.5, outlier.shape=NA) +
  geom_jitter(width=0.1, size=1, alpha=0.8) +
  scale_color_manual(values=c("grey80", "midnightblue", "darkred")) +
  scale_fill_manual(values=c("grey80", "darkred")) +
  scale_alpha_manual(values=c(0.001, 0.4)) +
  xlab("Line name") +
  ylab("Grain Yield (Kg/ha)") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  facet_wrap(~year)
# Save plot
ggsave("plot_level_spatially_adjusted_by_year.png", width=22, dpi=300)

### Plot level by location
ggplot(data=df.pl, aes(x=reorder(line_name, yield_kgha_use, FUN=median), y=as.numeric(yield_kgha_use), color=linetype, fill=seq_status, alpha=seq_status)) +
  geom_boxplot(width=0.5, outlier.shape=NA) +
  geom_jitter(width=0.1, size=1, alpha=0.8) +
  scale_color_manual(values=c("grey80", "midnightblue", "darkred")) +
  scale_fill_manual(values=c("grey80", "darkred")) +
  scale_alpha_manual(values=c(0.001, 0.4)) +
  xlab("Line name") +
  ylab("Grain Yield (Kg/ha)") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  facet_wrap(~location_name)
# Save plot
ggsave("plot_level_spatially_adjusted_by_location.png", width=22, dpi=300)
