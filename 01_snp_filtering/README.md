# SNP filtering

This directory contains all scripts used to prepare and filter VCF files.

## Dependencies

| Dependency | Version | Use |
| ---------- | ------- | --- |
| bcftools | 1.9 | checks, filtering, and indexing VCFs |

## Methods

When we got the 10x Genomics VCFs output from `longranger`, three of the samples had incorrect sample names that needed to be corrected. The only sample that did not need to be fixed was morex-sample2. Before working with the VCFs, I used bcftools reheader to change the sample names within the VCF files.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering
for i in reheader*.job; do qsub $i; done
```

**Step 0:** Prepare the VCFs

Do some checks to see which files have REF allele mismatches. We repeat the following check for all samples:

```bash
# In dir: ~/Shared/Projects/Mutant_Barley/longranger/M01-3-3-self5
REF=/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules_parts.fasta
# dels
bcftools norm -c w -f ${REF} -o temp_ref_check_M01-3-3_dels.vcf -O v M01-3-3_dels.vcf.gz >& temp_ref_check_warn_M01-3-3_dels.log
# Check unique warnings
cut -f 1 temp_ref_check_warn_M01-3-3_dels.log | sort -u
Lines   total/split/realigned/skipped:
NON_ACGTN_ALT

# Large SVs
bcftools norm -c w -f ${REF} -o temp_ref_check_M01-3-3_large_svs.vcf -O v M01-3-3_large_svs.vcf.gz >& temp_ref_check_warn_M01-3-3_large_svs.log
# Check unique warnings
cut -f 1 temp_ref_check_warn_M01-3-3_large_svs.log | sort -u
Lines   total/split/realigned/skipped:
NON_ACGTN_ALT
REF_MISMATCH
# Check how many REF_MISMATCH
grep "REF_MISMATCH" temp_ref_check_warn_M01-3-3_large_svs.log | wc -l
5926

# Phased variants
bcftools norm -c w -f ${REF} -o temp_ref_check_M01-3-3_phased_variants.vcf -O v M01-3-3_phased_variants.vcf.gz >& temp_ref_check_warn_M01-3-3_phased_variants.log
# Check unique warnings
cut -f 1 temp_ref_check_warn_M01-3-3_phased_variants.log | sort -u
Lines   total/split/realigned/skipped:

# Cleanup unnecessary files
rm temp_ref_check_*.vcf
```

After the checks, we found that the dels and phased variants vcf files do not have any REF mismatches (at the same position for multiple samples, the REF allele differed), only the large SVs vcf has REF mismatches. This is a bug in 10x Genomics longranger v2.2.2. We will exclude positions that have REF mismatches. NOTE: these largeSVs vcfs cannot be corrected with `bcftools norm -c wxs` because of the non-standard format used to report ALT alleles for deletions and large structural variants. If we used `bcftools norm -c wxs`, many of the positions would be removed due to the non-standard format used.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering
# Mutated lines
qsub -t 1-3 fix_ref_allele_large_svs.job
# morex-sample2
qsub -t 1 fix_ref_allele_large_svs-morex.job
```

For the phased variants vcf files relative to **Morex v2**, there were no REF allele mismatches. We also did a check to see if there were any missing REF alleles and we found none, so we can skip running the script below called `fix_ref_allele_phased_variants.job` since it only applied when we mapped to **Morex v1**.

```bash
# Check for missing REF alleles
bcftools filter -i 'REF == "."' M01-3-3_phased_variants.vcf.gz | less -S
bcftools filter -i 'REF == "."' M20-2-2_phased_variants.vcf.gz | less -S
bcftools filter -i 'REF == "."' M29-2-2_phased_variants.vcf.gz | less -S
bcftools filter -i 'REF == "."' morex-sample2_phased_variants.vcf.gz | less -S
```

(SKIP for Morex v2) For the phased variants vcf files relative to **Morex v1**, there were no REF allele mismatches, but there were many positions with missing REF alleles (this is also a bug in 10x Genomics longranger v2.2.2). We will fix any positions with missing REF allele.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering
# IMPORTANT: this step was only necessary for when mapped to Morex v1
qsub -t 1-4 fix_ref_allele_phased_variants.job
```

We want to focus on single bp changes and will remove any complex variants. We are doing this because it is difficult to say much about complex variants especially when we can't "place" them. Also, we ran into errors when using bcftools to merge the mutated morex vcf files. The complex variants seemed to be causing the issue.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering
# Submit a job to keep a record of complex variants and filter out complex variants
# Mutated lines
qsub -t 1-3 remove_complex_variants.job
# morex-sample2
qsub -t 1 remove_complex_variants-morex.job
```

Concatenate the following VCF files output from longranger:

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering
# Mutated lines
for i in concat_vcfs_M*.job; do qsub $i; done
# morex-sample2
qsub concat_vcfs_morex-sample2.job
```

**Step 1:** Merge together mutated line VCF files for 3 concatenated samples output from Longranger (keep Morex Longranger output separate) using bcftools.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering
qsub merge_vcfs.job

# Sort merged VCF files
qsub sort_vcfs.job
```

**Step 2:** Pull mutated and Morex VCFs into IGV to see how much heterogeneity there is

**Step 3:** Remove (subtract) everything that intersects between morex and mutated lines (using Longranger’s filtering)

Identify differences between 10x Morex and Morex references (outputs a BED file containing differences from reference).

```bash
~/GitHub/Barley_Mutated/01_snp_filtering/variant_filtering-morex-sample2.sh
```

Filter mutated lines by quality metrics (using 10x Genomics custom filters and DP), exclude sitest that differ between 10x Morex and Morex reference, and exclude non-unique variants (i.e., variants that are present in more than 1 of the mutated lines).

```bash
~/GitHub/Barley_Mutated/01_snp_filtering/variant_filtering-mut_3_lines.sh
```
