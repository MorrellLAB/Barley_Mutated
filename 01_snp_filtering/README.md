# SNP filtering

This directory contains all scripts used to prepare and filter VCF files.

## Dependencies

| Dependency | Version | Use |
| ---------- | ------- | --- |
| bcftools | 1.9 | checks, filtering, and indexing VCFs |

## Methods

When we got the 10x Genomics VCFs output from `longranger`, three of the samples had incorrect sample names that needed to be corrected. The only sample that did not need to be fixed was morex-sample2. Before working with the VCFs, I used bcftools reheader to change the sample names within the VCF files.

```bash
# In dir: ~/GitHub/Barley_Mutated/sequence_processing
for i in reheader*.job; do qsub $i; done
```

**Step 0:** Prepare the VCFs

Do some checks to see which files have REF allele mismatches. We repeat the following check for all samples:

```bash
# In dir: ~/Shared/Projects/Mutant_Barley/longranger/M01-3-3-self5
REF=/home/morrellp/liux1299/Shared/References/Reference_Sequences/Barley/Morex/barley_RefSeq_v1.0/barley_pseudomolecules_parts.fa
# dels
bcftools norm -c w -f ${REF} -o temp_ref_check_M01-3-3_dels.vcf -O v M01-3-3_new_samp_name_dels.vcf.gz >& temp_ref_check_warn_M01-3-3_dels.log
# Check unique warnings
cut -f 1 temp_ref_check_warn_M01-3-3_dels.log | sort -u
Lines   total/split/realigned/skipped:
NON_ACGTN_ALT
NON_ACGTN_REF

# Large SVs
bcftools norm -c w -f ${REF} -o temp_ref_check_M01-3-3_large_svs.vcf -O v M01-3-3_new_samp_name_large_svs.vcf.gz >& temp_ref_check_warn_M01-3-3_large_svs.log
# Check unique warnings
cut -f 1 temp_ref_check_warn_M01-3-3_large_svs.log | sort -u
Lines   total/split/realigned/skipped:
NON_ACGTN_ALT
NON_ACGTN_REF
REF_MISMATCH

# Phased variants
bcftools norm -c w -f ${REF} -o temp_ref_check_M01-3-3_phased_variants.vcf -O v M01-3-3_phased_variants.vcf.gz >& ref_check_warn_M01-3-3_phased_variants.log
# Check unique warnings
cut -f 1 ref_check_warn_M01-3-3_phased_variants.log | sort -u
Lines   total/split/realigned/skipped:
NON_ACGTN_REF

# Cleanup unnecessary files
rm temp_ref_check_*.vcf
```

After the checks, we found that the dels and phased variants vcf files do not have any REF mismatches (at the same position for multiple samples, the REF allele differed), only the large SVs vcf has REF mismatches. This is a bug in 10x Genomics longranger v2.2.2. We will exclude positions that have REF mismatches. NOTE: these cannot be corrected with `bcftools norm -c wxs` because of the non-standard format used to report ALT alleles for deletions and large structural variants. If we used `bcftools norm -c wxs`, many of the positions would be removed due to the non-standard format used.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering
qsub -t 1-4 fix_ref_allele_large_svs.job
```

For the phased variants vcf files, there were no REF allele mismatches, but there were many positions with missing REF alleles (this is also a bug in 10x Genomics longranger v2.2.2). We will fix any positions with missing REF allele.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_sequence_processing
qsub -t 1-4 fix_ref_allele_phased_variants.job
```

We want to focus on single bp changes and will remove any complex variants. We are doing this because it is difficult to say much about complex variants especially when we can't "place" them. Also, we ran into errors when using bcftools to merge the mutated morex vcf files. The complex variants seemed to be causing the issue.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_sequence_processing
# Submit a job to keep a record of complex variants and filter out complex variants
qsub -t 1-4 remove_complex_variants.job
```

Concatenate the following VCF files output from longranger:

```bash
# In dir: ~/GitHub/Barley_Mutated/sequence_processing
for i in concat*.job; do qsub $i; done
```

**Step 1:** Merge together mutated line VCF files for 3 concatenated samples output from Longranger (keep Morex Longranger output separate) using bcftools.

```bash
# In dir: ~/GitHub/Barley_Mutated/sequence_processing
qsub merge_vcfs.job
```

**Step 2:** Pull mutated and Morex VCFs into IGV to see how much heterogeneity there is

**Step 3:** Remove (subtract) everything that intersects between morex and mutated lines (using Longranger’s filtering)
