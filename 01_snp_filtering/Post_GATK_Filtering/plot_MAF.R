fp1 <- "~/Downloads/temp_msi/mut8_and_hybrid_barley_snps_poly_and_biallelic.maf"
fp2 <- "~/Downloads/temp_msi/mut8_and_hybrid_barley_snps_biallelic.maf"
fp3 <- "~/Downloads/temp_msi/mut8_and_hybrid_barley_snps_biallelic.noRepeatOverlap.noRefNs.maf"

df <- read.delim(fp1, header = TRUE, sep = '\t')
df_filt1 <- read.delim(fp2, header = TRUE, sep = '\t')
df_filt2 <- read.delim(fp3, header = TRUE, sep = '\t')


# Test plot to experiment with bin size
hist(df$MAF, breaks = 100, xlab="MAF", main="Mut8 and Hybrid poly and bi-allelic Morex V3 (100 bins)")
# Enable multi-panel plots
par(mfrow=c(2,2), mar=c(4.1, 4.1, 4.1, 2.1))
hist(df$MAF, breaks = 60, xlab="MAF", main="Mut8 and Hybrid polymorphic SNPs Morex V3")

# Filtered Morex v3
hist(df_filt1$MAF, breaks = 60, xlab="MAF", main="Mut8 and Hybrid filtered 1 Morex V3")
hist(df_filt2$MAF, breaks = 60, xlab="MAF", main="Mut8 and Hybrid filtered 2 Morex V3")
