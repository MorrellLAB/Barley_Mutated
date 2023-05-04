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
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        title=element_text(size=20),
        axis.text.x=element_text(angle=45, hjust=1),
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
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        title=element_text(size=20),
        axis.text.x=element_text(angle=45, hjust=1),
        plot.title=element_text(hjust=0.5))

grid.arrange(p1.sa.avg_yield, p2.raw.avg_yield, nrow=2)
combo.avg_yield <- arrangeGrob(p1.sa.avg_yield, p2.raw.avg_yield, nrow=2)
# Save plot
ggsave("trial_avg_spatially_adjusted_vs_raw.jpg", combo.avg_yield, width=13, height=10, units="in", dpi=300)

#-----------------
# Plot level average
df.pl <- read.csv(file=fp_plot_level_avg, header=TRUE, skip=1)
# Remove "FILLER" rows
df.pl <- df.pl %>% filter(!grepl("FILLER", line_name))

# Convert yield to numeric
df.pl$yield_bua_use <- as.numeric(df.pl$yield_bua_use)
df.pl$yield_kgha_use <- as.numeric(df.pl$yield_kgha_use)
df.pl$yield_kgha <- as.numeric(df.pl$yield_kgha)
# Convert other traits of interest to numeric
df.pl$heading_dap_use <- as.numeric(df.pl$heading_dap_use)
df.pl$height_cm_use <- as.numeric(df.pl$height_cm_use)

# Highlight mutated lines
df.pl <- df.pl %>% mutate(linetype=case_when((line_name == "MOREX_W2017") ~ "Parent",
                                         (line_name %in% check_lines) ~ "Check",
                                         !(line_name %in% check_lines) ~ "Sodium azide"))

# Shade sequenced mutated lines
df.pl <- df.pl %>% mutate(seq_status=ifelse(line_name %in% seq_mut_lines, "Sequenced", "Not Sequenced"))

# Get number of replicates in each check line
df.pl %>%
  filter(linetype == "Check") %>%
  select(line_name, rep, linetype) %>%
  group_by(line_name) %>%
  summarise(total_reps = sum(rep))

# Quick count number of mutated samples
df.pl %>% filter(linetype == "Sodium azide") %>% select(line_name) %>% arrange(line_name) %>% unique() %>% nrow()

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

# Create column of short names
uniq_line_names <- sort(unique(df.pl$line_name))
short_names_df <- data.frame(short_name=sort(vapply(strsplit(uniq_line_names, split="-"), '[', 1, FUN.VALUE=character(1))), line_name=uniq_line_names)

# Combine trait means into single df
# And add short names
df.pl.trait_means <- left_join(short_names_df, left_join(tmp.df.pl.yield, left_join(tmp.df.pl.hdap, tmp.df.pl.height, by="line_name"), by="line_name"), by="line_name")

############################
# Calculate average diminution in yield relative to Morex W2017 parent
morex_parent_yield <- df.pl %>%
  select(c(line_name, yield_kgha_use, yield_bua_use, yield_kgha)) %>%
  filter(line_name == "MOREX_W2017") %>%
  group_by(line_name) %>%
  summarise(across(everything(), mean))

morex_parent_hdap <- df.pl %>%
  select(c(line_name, heading_dap_use)) %>%
  filter(line_name == "MOREX_W2017") %>%
  filter(complete.cases(.)) %>%
  group_by(line_name) %>%
  summarise(across(everything(), mean))

morex_parent_height <- df.pl %>%
  select(c(line_name, height_cm_use)) %>%
  filter(line_name == "MOREX_W2017") %>%
  filter(complete.cases(.)) %>%
  group_by(line_name) %>%
  summarise(across(everything(), mean))

# Combine trait means into single df
df.morex.trait_means <- left_join(morex_parent_yield, left_join(morex_parent_hdap, morex_parent_height, by="line_name"), by="line_name")

df.pl.trait_means <- df.pl.trait_means %>%
  mutate(dim_yield_kgha_use = yield_kgha_use - df.morex.trait_means$yield_kgha_use,
         dim_heading_dap_use = mean_heading_dap_use - df.morex.trait_means$heading_dap_use,
         dim_height_cm_use = mean_height_cm_use - df.morex.trait_means$height_cm_use)

# Highlight treated vs untreated lines
df.pl.trait_means <- df.pl.trait_means %>% mutate(linetype=case_when((line_name == "MOREX_W2017") ~ "Untreated",
                                             (line_name %in% check_lines) ~ "Untreated",
                                             !(line_name %in% check_lines) ~ "Treated"))

# Add info on sequencing status of mutated lines
df.pl.trait_means <- df.pl.trait_means %>% mutate(seq_status=ifelse(line_name %in% seq_mut_lines, "Sequenced", "Not Sequenced"))

# Save df to file
df.pl.trait_means %>% write_csv(file="diminution_in_trait_means_per_sample.csv")

# Calculate average diminution
df.pl.trait_means %>%
  filter(linetype=="Treated") %>%
  summarise(mean_yield_kgha_use=mean(yield_kgha_use), mean_dim_yield_kgha_use=mean(dim_yield_kgha_use))

df.pl.trait_means %>%
  group_by(linetype) %>%
  summarize(cor_linetype=cor(dim_yield_kgha_use, dim_height_cm_use))

############################
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
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        title=element_text(size=20),
        axis.text.x=element_text(angle=45, hjust=1),
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
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        title=element_text(size=20),
        axis.text.x=element_text(angle=45, hjust=1),
        plot.title=element_text(hjust=0.5))

grid.arrange(p1.sa.yield, p2.raw.yield, nrow=2)
combo.plot_yield <- arrangeGrob(p1.sa.yield, p2.raw.yield, nrow=2)
# Save plot
ggsave("plot_level_spatially_adjusted_vs_raw.jpg", combo.plot_yield, width=13, height=10, units="in", dpi=300)

# Spatially adjusted only
ggplot(data=df.pl, aes(x=reorder(line_name, yield_kgha_use, FUN=median), y=as.numeric(yield_kgha_use), color=linetype, fill=seq_status, alpha=seq_status)) +
  geom_boxplot(width=0.6, outlier.shape=NA) +
  geom_jitter(width=0.2, size=1, alpha=0.8) +
  scale_color_manual(name="Sample Type", values=c("grey80", "midnightblue", "darkred")) +
  scale_fill_manual(name="Sequenced Mutated Lines", values=c("grey80", "darkred")) +
  scale_alpha_manual(name="Sequenced Mutated Lines", values=c(0.001, 0.4)) +
  xlab("Line name") +
  ylab("Grain Yield (Kg/ha)") +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        axis.text.x=element_text(angle=45, hjust=1),
        plot.title=element_text(hjust=0.5))
# Save just spatially adjusted
ggsave("plot_level_spatially_adjusted.jpg", width=14, height=8, dpi=300)

### Plot level by year
ggplot(data=df.pl, aes(x=reorder(line_name, yield_kgha_use, FUN=median), y=as.numeric(yield_kgha_use), color=linetype, fill=seq_status, alpha=seq_status)) +
  geom_boxplot(width=0.6, outlier.shape=NA) +
  geom_jitter(width=0.1, size=1, alpha=0.8) +
  scale_color_manual(name="Sample Type", values=c("grey80", "midnightblue", "darkred")) +
  scale_fill_manual(name="Sequenced Mutated Lines", values=c("grey80", "darkred")) +
  scale_alpha_manual(name="Sequenced Mutated Lines", values=c(0.001, 0.4)) +
  xlab("Line name") +
  ylab("Grain Yield (Kg/ha)") +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=32),
        legend.text=element_text(size=28),
        legend.title=element_text(size=28, face="bold"),
        strip.text=element_text(size=36, face="bold"),
        legend.position="bottom", legend.box="vertical",
        axis.text.x=element_text(angle=45, hjust=1)) +
  facet_wrap(~year)
# Save plot
ggsave("plot_level_spatially_adjusted_by_year.jpg", width=27, height=12, units="in", dpi=300)

### Plot level by location
ggplot(data=df.pl, aes(x=reorder(line_name, yield_kgha_use, FUN=median), y=as.numeric(yield_kgha_use), color=linetype, fill=seq_status, alpha=seq_status)) +
  geom_boxplot(width=0.6, outlier.shape=NA) +
  geom_jitter(width=0.1, size=1, alpha=0.8) +
  scale_color_manual(name="Sample Type", values=c("grey80", "midnightblue", "darkred")) +
  scale_fill_manual(name="Sequenced Mutated Lines", values=c("grey80", "darkred")) +
  scale_alpha_manual(name="Sequenced Mutated Lines", values=c(0.001, 0.4)) +
  xlab("Line name") +
  ylab("Grain Yield (Kg/ha)") +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=32),
        legend.text=element_text(size=28),
        legend.title=element_text(size=28, face="bold"),
        strip.text=element_text(size=36, face="bold"),
        legend.position="bottom", legend.box="vertical",
        axis.text.x=element_text(angle=45, hjust=1)) +
  facet_wrap(~location_name)
# Save plot
ggsave("plot_level_spatially_adjusted_by_location.jpg", width=27, height=12, units="in", dpi=300)
