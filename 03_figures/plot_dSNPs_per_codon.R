#!/usr/bin/env Rscript

library(tidyr)
library(ggplot2)
library(gridExtra)

# Generate plots of dSNPs per codon
# Inspired from Tom Kono's dSNP per codon plotting script:
# https://github.com/MorrellLAB/Deleterious_GP/blob/master/Analysis_Scripts/Plotting/Plot_Deleterious_per_Codon.R

# User provided input arguments
out_dir_fp <- "~/Dropbox/Projects/Barley_Mutated/analyses/bad_mutations/plots"

# Centromere and pericentromere
centromere_fp <- "~/GitHub/morex_reference/morex_v3/pericentromere/MorexV3_centromere_positions.tsv"
pericentromere_fp <- "~/GitHub/morex_reference/morex_v3/pericentromere/pericentromere_physPos.txt"

# Codons per window from GFF
mRNA_codons_per_win_fp <- "~/Dropbox/Projects/Barley_Mutated/analyses/bad_mutations/visualization/mRNA_10Mbp_num_codons_per_window.txt"
cds_codons_per_win_fp <- "~/Dropbox/Projects/Barley_Mutated/analyses/bad_mutations/visualization/CDS_10Mbp_num_codons_per_window.txt"

# Mutated
mut_del_snps_per_win_fp <- "~/Dropbox/Projects/Barley_Mutated/analyses/bad_mutations/visualization/mut_snps_private_num_del_snps_per_10Mbp_window.txt"
mut_tol_snps_per_win_fp <- "~/Dropbox/Projects/Barley_Mutated/analyses/bad_mutations/visualization/mut_snps_private_num_tol_snps_per_10Mbp_window.txt"

# Rare
rare_del_snps_per_win_fp <- "~/Dropbox/Projects/Barley_Mutated/analyses/bad_mutations/visualization/hybrid13_rare_snps_num_del_snps_per_10Mbp_window.txt"
rare_tol_snps_per_win_fp <- "~/Dropbox/Projects/Barley_Mutated/analyses/bad_mutations/visualization/hybrid13_rare_snps_num_tol_snps_per_10Mbp_window.txt"

# Common
comm_del_snps_per_win_fp <- "~/Dropbox/Projects/Barley_Mutated/analyses/bad_mutations/visualization/hybrid13_common_snps_num_del_snps_per_10Mbp_window.txt"
comm_tol_snps_per_win_fp <- "~/Dropbox/Projects/Barley_Mutated/analyses/bad_mutations/visualization/hybrid13_common_snps_num_tol_snps_per_10Mbp_window.txt"

#------------
# Set working directory for outputting plots
setwd(out_dir_fp)

# Read in data
read_num_per_win <- function(fp) {
  df <- read.delim(fp, header=TRUE, sep="\t")
  # Cleanup column name
  colnames(df)[1] = "chr"
  return(df)
}

prep_codons_per_win <- function(df_codons_per_win, df_del_vs_tol, pcent_df) {
  # Combine del vs tol with per codon counts
  df_merged <- merge(df_codons_per_win, df_del_vs_tol, by=c("chr", "start_pos", "end_pos"))
  # Prepare window midpoint positions
  df_merged$midp <- (df_merged$start_pos + df_merged$end_pos)/2
  # Convert from bp to Mbp
  df_merged$midp_mbp <- df_merged$midp/1000000
  # Calculate dSNP per codon and tolSNP per codon
  df_merged$snp_per_codon <- df_merged$num_snp/df_merged$num_codons
  # Add pericentromere info to df
  df_merged <- merge(df_merged, pcent_df, by="chr", all.x=TRUE)
  return(df_merged)
}

prep_data_format <- function(del_snps_per_win_fp, tol_snps_per_win_fp, df_codons_per_win, pericentromere_df) {
  df_dSNPs_per_win <- read_num_per_win(del_snps_per_win_fp)
  df_tolSNPs_per_win <- read_num_per_win(tol_snps_per_win_fp)
  
  # Rename columns for concatenating
  colnames(df_dSNPs_per_win)[5] = "num_snp"
  colnames(df_tolSNPs_per_win)[5] = "num_snp"
  
  # Combine del and tol df
  df_del_vs_tol <- rbind(df_dSNPs_per_win, df_tolSNPs_per_win)
  # Sort dataframe
  df_del_vs_tol <- df_del_vs_tol[order(df_del_vs_tol$chr, df_del_vs_tol$start_pos), ]
  
  # Codons per win
  df_codons <- prep_codons_per_win(df_codons_per_win, df_del_vs_tol, pericentromere_df)
  return(df_codons)
}

# Read in centromere and pericentromere positions
centromere_df <- read.delim(centromere_fp, header=FALSE, sep=" ")
colnames(centromere_df) <- c("chr", "cent_pos")
# Convert bp to Mbp
centromere_df$cent_pos_mbp <- centromere_df$cent_pos/1000000

pericentromere_df <- read.delim(pericentromere_fp, header=FALSE, sep="\t")
colnames(pericentromere_df) <- c("chr", "pcent_start", "pcent_end")
# Convert bp to Mbp
pericentromere_df$pcent_start_mbp <- pericentromere_df$pcent_start/1000000
pericentromere_df$pcent_end_mbp <- pericentromere_df$pcent_end/1000000

# Read in data
# Codons per window
df_mRNA_codons_per_win <- read_num_per_win(mRNA_codons_per_win_fp)
df_cds_codons_per_win <- read_num_per_win(cds_codons_per_win_fp)

# Define colors for Deleterious vs Tolerated
del_col <- '#d7191c'
tol_col <- '#fdae61'

###############
### Mutated ###
###############
df_mut_mRNA <- prep_data_format(mut_del_snps_per_win_fp, mut_tol_snps_per_win_fp, df_mRNA_codons_per_win, pericentromere_df)
df_mut_cds <- prep_data_format(mut_del_snps_per_win_fp, mut_tol_snps_per_win_fp, df_cds_codons_per_win, pericentromere_df)

# Plot dSNPs per codon: Deleterious vs Tolerated
# mRNA
ggplot(df_mut_mRNA, aes(x=midp_mbp, y=snp_per_codon, color=factor(del_vs_tol))) +
  #geom_rect(aes(xmin=pcent_start_mbp, xmax=pcent_end_mbp, ymin=-Inf, ymax=Inf), alpha=0.1, fill="grey", color="grey", na.rm=TRUE) +
  geom_vline(data=centromere_df, aes(xintercept=cent_pos_mbp), colour="darkgrey", size=1) +
  geom_point(size=3, alpha=0.8) +
  theme_bw() +
  xlab("Position (Mbp)") +
  ylab("SNPs per Codon") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=16),
        panel.grid=element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(size=16, face="bold")) +
  facet_wrap(~chr, ncol=1, scales='free') +
  scale_color_manual(values=c(del_col, tol_col), name="") +
  scale_x_continuous(limits=c(0, 700), n.breaks=24) +
  scale_y_continuous(limits=c(0, max(df_mut_mRNA$snp_per_codon)+(max(df_mut_mRNA$snp_per_codon)*0.2)))
# Save plot
ggsave(filename="dSNPs_per_codon_mRNA-mut.jpg", width=13, height=12, units="in", dpi=300)

# Plot dSNPs per codon: Deleterious vs Tolerated
# CDS
ggplot(df_mut_cds, aes(x=midp_mbp, y=snp_per_codon, color=factor(del_vs_tol))) +
  geom_vline(data=centromere_df, aes(xintercept=cent_pos_mbp), colour="darkgrey", size=1) +
  geom_point(size=3, alpha=0.75) +
  theme_bw() +
  xlab("Position (Mbp)") +
  ylab("SNPs per Codon") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=16),
        panel.grid=element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(size=16, face="bold")) +
  facet_wrap(~chr, ncol=1, scales='free') +
  scale_color_manual(values=c(del_col, tol_col), name="") +
  scale_x_continuous(limits=c(0, 700), n.breaks=24) +
  scale_y_continuous(limits=c(0, max(df_mut_cds$snp_per_codon)+(max(df_mut_cds$snp_per_codon)*0.2)))
# Save plot
ggsave(filename="dSNPs_per_codon_cds-mut.jpg", width=13, height=12, units="in", dpi=300)

############
### Rare ###
############
df_rare_mRNA <- prep_data_format(rare_del_snps_per_win_fp, rare_tol_snps_per_win_fp, df_mRNA_codons_per_win, pericentromere_df)
df_rare_cds <- prep_data_format(rare_del_snps_per_win_fp, rare_tol_snps_per_win_fp, df_cds_codons_per_win, pericentromere_df)

# Plot dSNPs per codon: Deleterious vs Tolerated
# mRNA
ggplot(df_rare_mRNA, aes(x=midp_mbp, y=snp_per_codon, color=factor(del_vs_tol))) +
  geom_vline(data=centromere_df, aes(xintercept=cent_pos_mbp), colour="darkgrey", size=1) +
  geom_point(size=3, alpha=0.75) +
  theme_bw() +
  xlab("Position (Mbp)") +
  ylab("SNPs per Codon") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=16),
        panel.grid=element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(size=16, face="bold")) +
  facet_wrap(~chr, ncol=1, scales='free') +
  scale_color_manual(values=c(del_col, tol_col), name="") +
  scale_x_continuous(limits=c(0, 700), n.breaks=24) +
  scale_y_continuous(limits=c(0, max(df_rare_mRNA$snp_per_codon)+(max(df_rare_mRNA$snp_per_codon)*0.2)))
# Save plot
ggsave(filename="dSNPs_per_codon_mRNA-hybrid_rare.jpg", width=13, height=12, units="in", dpi=300)

# Plot dSNPs per codon: Deleterious vs Tolerated
# CDS
ggplot(df_rare_cds, aes(x=midp_mbp, y=snp_per_codon, color=factor(del_vs_tol))) +
  geom_vline(data=centromere_df, aes(xintercept=cent_pos_mbp), colour="darkgrey", size=1) +
  geom_point(size=3, alpha=0.75) +
  theme_bw() +
  xlab("Position (Mbp)") +
  ylab("SNPs per Codon") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=16),
        panel.grid=element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(size=16, face="bold")) +
  facet_wrap(~chr, ncol=1, scales='free') +
  scale_color_manual(values=c(del_col, tol_col), name="") +
  scale_x_continuous(limits=c(0, 700), n.breaks=24) +
  #scale_y_continuous(limits=c(0, max(df_rare_cds$snp_per_codon)+(max(df_rare_cds$snp_per_codon)*0.2)))
  scale_y_continuous(limits=c(0, 0.012))
# Save plot
ggsave(filename="dSNPs_per_codon_cds-hybrid_rare.jpg", width=13, height=12, units="in", dpi=300)

##############
### Common ###
##############
df_common_mRNA <- prep_data_format(comm_del_snps_per_win_fp, comm_tol_snps_per_win_fp, df_mRNA_codons_per_win, pericentromere_df)
df_common_cds <- prep_data_format(comm_del_snps_per_win_fp, comm_tol_snps_per_win_fp, df_cds_codons_per_win, pericentromere_df)

# Plot dSNPs per codon: Deleterious vs Tolerated
# mRNA
ggplot(df_common_mRNA, aes(x=midp_mbp, y=snp_per_codon, color=factor(del_vs_tol))) +
  geom_vline(data=centromere_df, aes(xintercept=cent_pos_mbp), colour="darkgrey", size=1) +
  geom_point(size=3, alpha=0.75) +
  theme_bw() +
  xlab("Position (Mbp)") +
  ylab("SNPs per Codon") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=16),
        panel.grid=element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(size=16, face="bold")) +
  facet_wrap(~chr, ncol=1, scales='free') +
  scale_color_manual(values=c(del_col, tol_col), name="") +
  scale_x_continuous(limits=c(0, 700), n.breaks=24) +
  scale_y_continuous(limits=c(0, max(df_common_mRNA$snp_per_codon)+(max(df_common_mRNA$snp_per_codon)*0.2)))
# Save plot
ggsave(filename="dSNPs_per_codon_mRNA-hybrid_common.jpg", width=13, height=12, units="in", dpi=300)

# Plot dSNPs per codon: Deleterious vs Tolerated
# CDS
ggplot(df_common_cds, aes(x=midp_mbp, y=snp_per_codon, color=factor(del_vs_tol))) +
  geom_vline(data=centromere_df, aes(xintercept=cent_pos_mbp), colour="darkgrey", size=1) +
  geom_point(size=3, alpha=0.75) +
  theme_bw() +
  xlab("Position (Mbp)") +
  ylab("SNPs per Codon") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=16),
        panel.grid=element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(size=16, face="bold")) +
  facet_wrap(~chr, ncol=1, scales='free') +
  scale_color_manual(values=c(del_col, tol_col), name="") +
  scale_x_continuous(limits=c(0, 700), n.breaks=24) +
  #scale_y_continuous(limits=c(0, max(df_common_cds$snp_per_codon)+(max(df_common_cds$snp_per_codon)*0.2)))
  scale_y_continuous(limits=c(0, 0.012))
# Save plot
ggsave(filename="dSNPs_per_codon_cds-hybrid_common.jpg", width=13, height=12, units="in", dpi=300)

#########################################
### Representative chromosomes figure ###
#########################################
# Prepare new df of only chr1H and 2H and add dataset labels
chr1H_2H_mut_cds <- df_mut_cds[df_mut_cds$chr == "chr1H" | df_mut_cds$chr == "chr2H", ]
chr1H_2H_mut_cds$dataset <- "Mutated"
chr1H_2H_rare_cds <- df_rare_cds[df_rare_cds$chr == "chr1H" | df_rare_cds$chr == "chr2H", ]
chr1H_2H_rare_cds$dataset <- "Rare"
chr1H_2H_comm_cds <- df_common_cds[df_common_cds$chr == "chr1H" | df_common_cds$chr == "chr2H", ]
chr1H_2H_comm_cds$dataset <- "Common"
# Combine datasets
df_cds_chr1H_2H <- rbind(chr1H_2H_mut_cds, chr1H_2H_rare_cds, chr1H_2H_comm_cds)

# Plot chr1H and 2H of combined datasets
ggplot(df_cds_chr1H_2H, aes(x=midp_mbp, y=snp_per_codon, color=factor(del_vs_tol))) +
  geom_vline(data=centromere_df[centromere_df$chr == "chr1H" | centromere_df$chr == "chr2H", ], aes(xintercept=cent_pos_mbp), colour="darkgrey", size=1) +
  geom_point(size=3, alpha=0.75) +
  theme_bw() +
  xlab("Position (Mbp)") +
  ylab("SNPs per Codon") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text=element_text(size=18),
        legend.position="top",
        panel.grid=element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(size=18, face="bold")) +
  facet_grid(factor(dataset, levels=c('Mutated', 'Rare', 'Common')) ~ chr, scales='free') +
  scale_color_manual(values=c(del_col, tol_col), name="") +
  scale_x_continuous(limits=c(0, 700), n.breaks=24)
# Save plot
ggsave(filename="dSNPs_per_codon_cds-chr1H_2H.jpg", width=18, height=8, units="in", dpi=300)
