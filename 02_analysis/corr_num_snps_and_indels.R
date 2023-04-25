#!/usr/bin/env Rscript

library(tidyverse)
library(ggpubr)

# Correlation between number of SNPs and indels
# Also generate plot of counts per line

# User provided input arguments
num_snps_fp = "~/Dropbox/Projects/Barley_Mutated/analyses/counts_tables/mut_private_per_sample_counts.del_vs_tol.nonsyn_vs_syn.txt"
num_indels_fp = "~/Dropbox/Projects/Barley_Mutated/analyses/counts_tables/mut_private_indels_per_sample_counts.txt"
out_dir = "~/Dropbox/Projects/Barley_Mutated/analyses/counts_tables"
# Correlation plots output prefix
out_prefix = "corr_mut_private_per_sample-num_snps_vs_indels"
# Barplot output file prefix
bp_out_prefix = "mut_private_per_sample"

#------------
num_snps = read.delim(num_snps_fp, header=TRUE, sep="\t")
num_indels = read.delim(num_indels_fp, header=TRUE, sep="\t")

var_counts <- left_join(num_snps, num_indels, by="sample")

var_counts.no10x <- var_counts %>%
  filter(sample != 'M01-3-3' & sample != 'M20-2-2' & sample != 'M29-2-2')

# Test for normality counts
shapiro.test(var_counts$num_snps)
# p-value = 0.0002737
shapiro.test(var_counts$num_indels)

# Check QQ plots
ggqqplot(var_counts$num_snps, ylab="Number of SNPs")
ggqqplot(var_counts$num_indels, ylab="Number of indels")

# Correlation test
cor.test(var_counts$num_snps, var_counts$num_indels, method="spearman", alternative="two.sided")

jpeg(filename=paste0(out_dir, '/', out_prefix, '.jpg'), width=7, height=6, units='in', res=300)
ggscatter(var_counts, x="num_snps", y="num_indels",
          add="reg.line", conf.int=TRUE,
          cor.coef=TRUE, cor.method="spearman",
          xlab="Number of SNPs",
          ylab="Number of indels")
dev.off()

# Test for normality counts no 10x Genomics samples
shapiro.test((var_counts.no10x$num_snps))
shapiro.test((var_counts.no10x$num_indels))

# Check QQ plots
ggqqplot(var_counts.no10x$num_snps, ylab="Number of SNPs")
ggqqplot(var_counts.no10x$num_indels, ylab="Number of indels")

ggscatter(var_counts.no10x, x="num_snps", y="num_indels",
          add="reg.line", conf.int=TRUE,
          cor.coef=TRUE, cor.method="spearman")

####################
# Plots - use counts table created above
# Define colors for Deleterious vs Tolerated
del_col <- '#d7191c'
tol_col <- '#fdae61'
nonsyn_col <- '#ffb703'
#syn_col <- '#126782'
syn_col <- '#2c7da0'

# Convert table to long format
vc_del_tol.long <- var_counts %>% 
  select(sample, total_del, total_tol) %>%
  pivot_longer(cols=c("total_del", "total_tol"), names_to="del_vs_tol", values_to="count") %>%
  arrange(desc(count))
# Generate labels for stacked barplot by computing cumulative sums
# vc_del_tol.long <- vc_del_tol.long %>%
#   group_by(sample) %>%
#   arrange(sample, desc(count)) %>%
#   arrange(match(sample, c("M39-2-3_18-28")), count) %>%
#   mutate(cum_count = cumsum(count) - 0.5 * count)

# Stacked barplot: del vs tol
ggplot(vc_del_tol.long, aes(x=sample, y=count, fill=del_vs_tol)) +
  geom_bar(position="stack", stat="identity", width=0.75) +
  theme_bw() +
  scale_fill_manual(values=c(del_col, tol_col), labels=c("Deleterious", "Tolerated")) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        plot.title=element_text(hjust=0.5),
        legend.title=element_blank()) +
  xlab("Sample Name") +
  ylab("SNP Count")
# Save plot to file
ggsave(paste0(out_dir, '/', bp_out_prefix, '-del_vs_tol.jpg'), width=10, height=6, dpi=300)

# Stacked barplot: syn vs nonsyn
# Convert table to long format
vc_nonsyn_syn.long <- var_counts %>% 
  select(sample, total_syn, total_nonsyn) %>%
  pivot_longer(cols=c("total_syn", "total_nonsyn"), names_to="nonsyn_vs_syn", values_to="count") %>%
  arrange(desc(count))

ggplot(vc_nonsyn_syn.long, aes(x=sample, y=count, fill=nonsyn_vs_syn)) +
  geom_bar(position="stack", stat="identity", width=0.75) +
  theme_bw() +
  scale_fill_manual(values=c(nonsyn_col, syn_col), labels=c("Nonsynonymous", "Synonymous")) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        plot.title=element_text(hjust=0.5),
        legend.title=element_blank()) +
  xlab("Sample Name") +
  ylab("SNP Count")
# Save Plot
ggsave(paste0(out_dir, '/', bp_out_prefix, '-nonsyn_vs_syn.jpg'), width=10, height=6, dpi=300)
 