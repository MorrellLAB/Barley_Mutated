# SNP filtering

This directory contains all scripts used to prepare and filter VCF files.

## Dependencies

| Dependency | Version | Use |
| ---------- | ------- | --- |
| bcftools | 1.9 | checks, filtering, and indexing VCFs |

## Methods

**Step 0:** Prepare the VCFs and uncallable regions

**VCFs**

Prepare large SVs vcf so low confidence large SV calls are excluded.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering
./prep_10x_large_svs_files.sh
```

Do some checks to see which files have REF allele mismatches (this is a known problem we have run into with output from 10x Genomics longranger outputs). We repeat the following check for all samples:

```bash
# Dependencies
module load bcftools/1.10.2
export PATH=${PATH}:/panfs/jay/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/01_snp_filtering

# Shared variables
REF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta"

# Morex-sample2
# In dir: ~/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/outs
check_10x_ref_allele_mismatches.sh \
    ${REF} \
    morex-sample2_dels.vcf.gz \
    morex-sample2_large_svs.calls.nochrUn.noDUP-UNK.vcf \
    morex-sample2_phased_variants.vcf.gz
# Output messages:
# Unique warnings for DELs:
# Lines   total/split/realigned/skipped:
# NON_ACGTN_ALT

# Unique warnings for Large SVs:
# Lines   total/split/realigned/skipped:
# NON_ACGTN_ALT

# Unique warnings for Phased Variants:
# Lines   total/split/realigned/skipped:

# Check how many REF_MISMATCH
grep "REF_MISMATCH" temp_ref_check_warn_morex-sample2_large_svs.log | wc -l
#    0

# M01
# In dir: ~/Projects/Mutant_Barley/longranger_morex_v3/M01-3-3/outs
check_10x_ref_allele_mismatches.sh \
    ${REF} \
    M01-3-3_dels.vcf.gz \
    M01-3-3_large_svs.calls.nochrUn.noDUP-UNK.vcf \
    M01-3-3_phased_variants.vcf.gz
# Check how many REF_MISMATCH
grep "REF_MISMATCH" temp_ref_check_warn_M01-3-3_large_svs.log | wc -l
#    2

# M20
# In dir: ~/Projects/Mutant_Barley/longranger_morex_v3/M20-2-2/outs
check_10x_ref_allele_mismatches.sh \
    ${REF} \
    M20-2-2_dels.vcf.gz \
    M20-2-2_large_svs.calls.nochrUn.noDUP-UNK.vcf \
    M20-2-2_phased_variants.vcf.gz
# Check how many REF_MISMATCH
grep "REF_MISMATCH" temp_ref_check_warn_M20-2-2_large_svs.log | wc -l
#    2

# M29
# In dir: ~/Projects/Mutant_Barley/longranger_morex_v3/M29-2-2/outs
check_10x_ref_allele_mismatches.sh \
    ${REF} \
    M29-2-2_dels.vcf.gz \
    M29-2-2_large_svs.calls.nochrUn.noDUP-UNK.vcf \
    M29-2-2_phased_variants.vcf.gz
# Check how many REF_MISMATCH
grep "REF_MISMATCH" temp_ref_check_warn_M29-2-2_large_svs.log | wc -l
#   1

# Cleanup unnecessary files
rm temp_ref_check_*.vcf
```

After the checks, we found that the dels and phased variants vcf files do not have any REF mismatches (at the same position for multiple samples, the REF allele differed), only the large SVs vcf has REF mismatches. This is a bug in 10x Genomics longranger v2.2.2. We will exclude positions that have REF mismatches. NOTE: these largeSVs vcfs cannot be corrected with `bcftools norm -c wxs` because of the non-standard format used to report ALT alleles for deletions and large structural variants. If we used `bcftools norm -c wxs`, many of the positions would be removed due to the non-standard format used.

Remove ref allele mismatches for large SVs.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering
# mutated lines (morex-sample2 didn't have any ref mismatches)
sbatch --array=0-2 fix_ref_allele_large_svs.sh
```

A majority of the REF allele mismatches in the large SVs VCF files were `LOWQ` and would have been removed anyways, there were only a handful or less of sites where they were marked as `PASS`.

Check for missing REF alleles:

```bash
# Check for missing REF alleles
bcftools filter -i 'REF == "."' M01-3-3_phased_variants.vcf.gz
bcftools filter -i 'REF == "."' M20-2-2_phased_variants.vcf.gz
bcftools filter -i 'REF == "."' M29-2-2_phased_variants.vcf.gz
bcftools filter -i 'REF == "."' morex-sample2_phased_variants.vcf.gz
```

For the phased variants vcf files relative to **Morex v3**, there were no REF allele mismatches. There were also no missing REF alleles, so we can skip running the script below called `fix_ref_allele_phased_variants.job` (in `deprecated_scripts`) since it only applied when we mapped to **Morex v1**.

We want to focus on single bp changes and will remove any complex variants from the phased variants VCF. It is difficult to say much about complex variants especially when we can't "place" them. Also, we ran into errors when using bcftools to merge the mutated morex vcf files. The complex variants seemed to be causing the issue.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering
# Submit a job to keep a record of complex variants and filter out complex variants
# morex-sample2 and mutated lines
sbatch --array=0-3 remove_complex_variants.sh
```

**Uncallable regions**

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering
./prep_uncallable_regions.sh
```

**Step 1:** Filter VCFs.

We'll keep the `dels.vcf.gz`, `large_svs.vcf.gz`, and `phased_variants.vcf.gz` separate for the first part of our vcf filtering because each file has different types of annotations we can use for filtering so it is easier to keep them separate for now.

Filter morex-sample2 (10x Genomics VCF and Nanopore VCF), morex PacBio VCF, and smoove VCF. The scripts below were used and serve as a log of the filtering commands run.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering
# morex-sample2 10x Genomics VCF
# Did some filtering, then sent the filtered DELs to SV-plaudit and scored them.
#   Scripts and steps for SV-plaudit are in the subdirectory `Samplot-Morex`
vcf_filtering-morex_10x.sh
vcf_filtering-morex_ont.sh
vcf_filtering-morex_pacbio.sh
vcf_filtering-smoove.sh
```

Filter WGS sample GATK variant calls, scripts are in subdirectory `Post_GATK_Filtering`.

**Step 2:** Pull Morex and mutated VCFs into IGV to see how much heterogeneity there is

Morex VCFs include:
- 10x Genomics morex-sample2
- smoove called on 10x Genomics morex-sample2 data
- ONT morex-sample2
- ONT 85x Morex (100 seedlings from Mascher et al. 2021)
- PacBio Morex (100 seedlings from Mascher et al. 2021)

Mutated VCFs include:
- 10x Genomics three mutated lines (M01, M20, and M29)
- smoove called on 10x Genomics three mutated lines data

Tune VCF filtering if necessary.

**Step 3:** Remove (subtract) everything that intersects between morex and mutated lines (using Longranger’s filtering)

Identify differences between our morex sample (a.k.a morex-sample2) and the Morex reference. We want to output a BED file containing differences between Morex and Morex reference.

Filter mutated lines by quality metrics (using 10x Genomics custom filters and DP), exclude sites that differ between 10x Morex and Morex reference, and exclude non-unique variants (i.e., variants that are present in more than 1 of the mutated lines). The script below was used and serve as a log of the filtering commands run.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering
vcf_filtering-mut_3_lines.sh
```

Exclude differences between 10x Genomics morex-sample2 and Morex reference for 8 WGS mutated lines.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/Post_GATK_Filtering
./exclude_diffs_from_ref_mut8wgs.sh
```

Before finalizing filtering, check 10-15 variants in IGV to see if we need to go back and tune our filtering a little more before proceeding.

Remove variants shared between 10x Genomics Morex line and mutant lines. So, we want variants private to the mutant lines that are NOT in the 10x Genomics Morex line.

```bash
# In dir: ~/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated/Filtered
module load bcftools/1.9
bcftools isec -p /scratch.global/liux1299 mut_3_lines_filtered_singletons_only_annotated_DEL.vcf.gz
# 0000.vcf contains records private to mut_3_lines_filtered_singletons_only_annotated_DEL.vcf.gz
cp /scratch.global/liux1299/0000.vcf mut_3_lines_filtered_singletons_only_annotated_DEL_de_novo_sites.vcf
```

**Step 4:** Identify variants private to each mutated sample

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering
./vcf_filtering-de_novo_mut.sh
```

## File Locations

Filtered VCFs (quality control/minimize errors):


Differences between Morex samples and Morex reference:

