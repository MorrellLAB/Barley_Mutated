# morex-sample2 (10x Genomics data) diversity estimate

Goal is to identify windows where there is higher diversity than expected (likely due to mapping issues) between morex-sample2 and the morex reference genome. SNPs and indels called in these regions will be excluded from down stream analyses.

---

### Prepare invariant sites for pixy

Run Haplotype Caller to generate GVCF file (since this was not provided as part of the Longranger pipeline).

```bash
# In dir: ~/sequence_handling
# Parallelize across chr parts (otherwise run time takes way too long)
./sequence_handling Haplotype_Caller ~/GitHub/Barley_Mutated/01_snp_filtering/morex-sample2_diversity/Config_morex-sample2

# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/morex-sample2_diversity
# Run Genotype_GVCFs
# With special options to include invariant sites (for pixy)
sbatch --array=0-13 Genotype_GVCFs-morex-sample2.sh
```

Concatenate split VCF files.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/morex-sample2_diversity
sbatch concat_split_vcfs.sh
```

Pull out invariant sites only.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/morex-sample2_diversity
sbatch select_invariant_sites.sh
```

Select chr1H only to get a sense of annotation distributions.

```bash
# In dir: ~/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/Genotype_GVCFs
module load bcftools/1.10.2
module load datamash_ML/1.3

bcftools view --regions "chr1H_part1,chr1H_part2" morex-sample2.invariant_sites.vcf.gz -O z -o ~/scratch/morex-sample2_chr1H.invariant_sites.vcf.gz
zgrep -v "#" ~/scratch/morex-sample2_chr1H.invariant_sites.vcf.gz | cut -f 10 | cut -d':' -f 4 > ~/scratch/temp_RGQ.txt
# Downsample (too many variants takes a very long time)
grep -v "0/0" ~/scratch/temp_RGQ.txt | grep -v -e '^$' | awk 'NR % 5 == 0' > ~/scratch/temp_RGQ_downsampled.txt
# Get some summaries
datamash -H mean 1 q1 1 median 1 q3 1 iqr 1 sstdev 1 < ~/scratch/temp_RGQ_downsampled.txt
#mean(3)	q1(3)	median(3)	q3(3)	iqr(3)	sstdev(3)
#73.319554808213	52	81	100	48	33.579479542322
datamash -H perc:5 1 < ~/scratch/temp_RGQ_downsampled.txt
#perc:5(3)
#0
datamash -H perc:10 1 < ~/scratch/temp_RGQ_downsampled.txt
#perc:10(3)
#19
```

Filter invariant sites.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/morex-sample2_diversity
sbatch filter_invariant_sites.sh

# Check if thre are overlaps (discordant) with 10x Genomics Longranger SNPs/INDELs
sbatch intersect_vcf.sh
```

Remove indels form invariant sites (only need SNPs for pixy).

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/morex-sample2_diversity
sbatch remove_indels_invariant.sh

# Check number of indels
# In dir: ~/scratch/temp_morex-sample2
grep -v "#" morex-sample2.invariant_sites.filtDP-RGQ-QUAL.indels.vcf | wc -l
#789484
```

Remove discordant sites for the purposes of running pixy. In the filtered morex-sample2 invariant sites VCF file, there were a total of 558 snp or indel sites where the genotype call by Longranger differs from GATK (both VCFs are filtered). At these sites, GATK calls the site invariant but Longranger calls a variant at the site. These are a small enough number of sites that they will have a limited effect at the whole genome level, so we'll exclude them from pixy.

```bash
# In dir: ~/scratch/temp_morex-sample2
# Convert VCF to BED for sites to exclude
# Note: lines in the .txt files below are in VCF format, just missing header lines
module load bedops_ML/2.4.38
vcf2bed < temp_morex-sample2.10x_snp_x_invariant.txt | cut -f 1-3 > temp_morex-sample2.10x_snp_x_invariant.bed
vcf2bed < temp_morex-sample2.10x_indel_x_invariant.txt | cut -f 1-3 > temp_morex-sample2.10x_indel_x_invariant.bed
# Concatenate sites and sort
cat temp_morex-sample2.10x_snp_x_invariant.bed temp_morex-sample2.10x_indel_x_invariant.bed | sort -u -k1,1 -k2,2n > temp_morex-sample2.10x_snp_indel_x_invariant.bed

# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/morex-sample2_diversity
# Exclude discordant sites
sbatch exclude_discordant.sh
```

Concatenate variant and invariant sites.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/morex-sample2_diversity
# Invariant sites file is too large, so split by chromosome parts
# Combine invariant and variant sites for each chrom part and sort
sbatch --array=0-13 concat_variant_and_invariant.sh
# Concatenate split sorted chrom parts
sbatch concat_split_vcfs.sh
```

Prepare callable regions based on depth to estimate diversity. Purpose here is to identify windows where diversity is higher than expected and exclude those SNPs.

```bash
# In dir: ~/Projects/Mutant_Barley/longranger_morex_v3/mosdepth_coverage
zgrep -a -E 'CALLABLE|HIGH_COVERAGE' morex-sample2.quantized.bed.gz > morex-sample2.quantized.CALLABLE.HIGH_COVERAGE.bed
```

Diversity estimates for morex.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/morex-sample2_diversity
sbatch run_pixy.sh
```
