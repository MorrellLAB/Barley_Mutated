# Plot of del vs tol amino acid changes

library(tidyverse)
library(ggplot2)
library(forcats)

# Set working directory to output plots to
setwd("~/Dropbox/Projects/Barley_Mutated/analyses/bad_mutations/plots/AA_change_del_vs_tol")

# Files should be VeP results merged with BAD_Mutations deleterious vs tolerated annotations with pseudomolecular positions added
mut_fp <- "~/Dropbox/Projects/Barley_Mutated/analyses/bad_mutations/predictions/mut_snps_private_deleterious_vs_tolerated.with_pseudo_pos.txt"
hybrid_rare_fp <- "~/Dropbox/Projects/Barley_Mutated/analyses/bad_mutations/predictions/hybrid13_rare_snps_deleterious_vs_tolerated.with_pseudo_pos.txt"
hybrid_comm_fp <- "~/Dropbox/Projects/Barley_Mutated/analyses/bad_mutations/predictions/hybrid13_common_snps_deleterious_vs_tolerated.with_pseudo_pos.txt"

#-----------------
# Load data and add data group label
mut_df <- read.table(mut_fp, header=TRUE) %>% add_column(Data_Group="Mutated")
rare_df <- read.table(hybrid_rare_fp, header=TRUE) %>% add_column(Data_Group="Rare")
comm_df <- read.table(hybrid_comm_fp, header=TRUE) %>% add_column(Data_Group="Common")
# Concatenate dfs
all_df <- rbind(mut_df, rare_df, comm_df)

calc_aa_prop <- function(df, group_name) {
    aa <- df %>%
        group_by(Amino_acids) %>%
        summarise(n_aa = n()) %>%
        mutate(prop_aa = n_aa / sum(n_aa)) %>%
        as.data.frame() %>%
        add_column(Data_Group=group_name)
    return(aa)
}

# Proportion is grouped by deleterious vs tolerated
# Ex: Num A/T Deleterious / Total Deleterious
calc_aa_del_vs_tol_prop <- function(df, group_name) {
    aa_del_vs_tol <- df %>%
        group_by(Amino_acids, Del_vs_Tol) %>%
        summarise(n_aa_del_tol = n()) %>%
        ungroup() %>%
        as.data.frame() %>%
        group_by(Del_vs_Tol) %>%
        mutate(prop_aa_del_tol = n_aa_del_tol / sum(n_aa_del_tol)) %>%
        as.data.frame() %>%
        add_column(Data_Group=group_name)
    return(aa_del_vs_tol)
}

# Count and calculate proportions
mut_aa <- calc_aa_prop(mut_df, group_name="Mutated")
mut_aa_del_vs_tol <- calc_aa_del_vs_tol_prop(mut_df, group_name="Mutated")

# Hybrid Rare and Common
rare_aa <- calc_aa_prop(rare_df, group_name="Rare")
rare_aa_del_vs_tol <- calc_aa_del_vs_tol_prop(rare_df, group_name="Rare")

comm_aa <- calc_aa_prop(comm_df, group_name="Common")
comm_aa_del_vs_tol <- calc_aa_del_vs_tol_prop(comm_df, group_name="Common")

# Concatenate
all_aa_only <- rbind(mut_aa, rare_aa, comm_aa)
# Reorder groups for plotting
all_aa_only$Data_Group <- factor(all_aa_only$Data_Group, levels=c("Mutated", "Rare", "Common"))

all_aa_del_vs_tol <- rbind(mut_aa_del_vs_tol, rare_aa_del_vs_tol, comm_aa_del_vs_tol)
# Reorder groups for plotting
all_aa_del_vs_tol$Data_Group <- factor(all_aa_del_vs_tol$Data_Group, levels=c("Mutated", "Rare", "Common"))

#------------------
# Generate plots
# Define colors for Deleterious vs Tolerated
del_col <- '#d7191c'
tol_col <- '#fdae61'
# Colors for Mut vs Rare vs Common
mut_color <- '#d22525'
rare_color <- '#cee3f8'
comm_color <- '#6898d6'

### Frequency Plots ###
# Mut AA counts del vs tol
mut_aa_order <- c(mut_aa %>% arrange(desc(n_aa)) %>% pull(Amino_acids))
# Custom order swap based on experimenting with plots
# Swap order of A/V and P/S
mut_aa_order[3] <- "P/S"
mut_aa_order[4] <- "A/V"

mut_reordered_df <- mut_df
mut_reordered_df$Amino_acids <- factor(mut_reordered_df$Amino_acids, levels=c(mut_aa_order))

ggplot(mut_reordered_df, aes(x=fct_rev(Amino_acids), fill=fct_rev(Del_vs_Tol))) +
    geom_bar() +
    scale_fill_manual(values=c(tol_col, del_col)) +
    xlab("Amino Acid Change") +
    ylab("Count") +
    coord_flip() +
    theme_classic() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16),
          legend.text=element_text(size=14),
          plot.title=element_text(hjust=0.5),
          legend.title=element_blank())
# Save plot to file
ggsave("counts_mut_AA_change_del_vs_tol.jpg", width=8, height=9, dpi=300)

# Rare AA counts del vs tol
ggplot(rare_df, aes(x=fct_rev(fct_infreq(Amino_acids)), fill=fct_rev(Del_vs_Tol))) +
    geom_bar() +
    scale_fill_manual(values=c(tol_col, del_col)) +
    xlab("Amino Acid Change") +
    ylab("Count") +
    coord_flip() +
    theme_classic() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=16),
          legend.text=element_text(size=14),
          plot.title=element_text(hjust=0.5),
          legend.title=element_blank())
# Save plot to file
ggsave("counts_rare_AA_change_del_vs_tol.jpg", width=12, height=22, dpi=300)

# Common AA counts del vs tol
ggplot(comm_df, aes(x=fct_rev(fct_infreq(Amino_acids)), fill=fct_rev(Del_vs_Tol))) +
    geom_bar() +
    scale_fill_manual(values=c(tol_col, del_col)) +
    xlab("Amino Acid Change") +
    ylab("Count") +
    coord_flip() +
    theme_classic() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=16),
          legend.text=element_text(size=14),
          plot.title=element_text(hjust=0.5),
          legend.title=element_blank())
# Save plot to file
ggsave("counts_common_AA_change_del_vs_tol.jpg", width=12, height=22, dpi=300)

### Proportions Plots ###
# Stacked bar chart mut/rare/comm del vs tol proportions (plotted separately b/c deleterious prop is relative to sum of all deleterious)
all_aa_del_prop <- all_aa_del_vs_tol %>% filter(Del_vs_Tol == "Deleterious")
all_aa_tol_prop <- all_aa_del_vs_tol %>% filter(Del_vs_Tol == "Tolerated")
# Deleterious proportions
ggplot(all_aa_del_prop, aes(x=fct_rev(Amino_acids), y=prop_aa_del_tol)) +
    geom_bar(position="stack", stat="identity", fill=del_col, width=0.8) +
    xlab("Amino Acid Change") +
    ylab("Proportion") +
    coord_flip() +
    facet_wrap(~ Data_Group, ncol=3) +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=16),
          legend.text=element_text(size=14),
          plot.title=element_text(hjust=0.5),
          legend.title=element_blank(),
          strip.text=element_text(size=16))
# Save plot to file
ggsave("prop_AA_change_deleterious.jpg", width=14, height=22, dpi=300)

# Tolerated proportions
ggplot(all_aa_tol_prop, aes(x=fct_rev(Amino_acids), y=prop_aa_del_tol)) +
    geom_bar(position="stack", stat="identity", fill=tol_col, width=0.8) +
    xlab("Amino Acid Change") +
    ylab("Proportion") +
    coord_flip() +
    facet_wrap(~ Data_Group, ncol=3) +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=16),
          legend.text=element_text(size=14),
          plot.title=element_text(hjust=0.5),
          legend.title=element_blank(),
          strip.text=element_text(size=16))
# Save plot to file
ggsave("prop_AA_change_tolerated.jpg", width=14, height=22, dpi=300)
