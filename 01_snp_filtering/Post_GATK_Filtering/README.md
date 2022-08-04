# Post GATK Filtering for WGS Illumina-paired end reads

We followed the GATK best practices pipeline for SNP calling implemented in [`sequence_handling`](https://github.com/MorrellLAB/sequence_handling) all the way through the handler `Variant_Recalibrator`. The post GATK filtering steps are documented below.

### Methods: SNPs

Keep only polymorphic sites.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/Post_GATK_Filtering
sbatch Filter_VCF_Poly.sh
```

Next, filter on proportion heterozygous genotypes, min DP, proportion missing, quality, GQ, and filter by allelic balance. Also filter to biallelic sites, also pull out triallelic sites for exploratory purposes only.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/Post_GATK_Filtering
Filter_VCF_Het_AB_DP_Miss_Qual_GQ.sh
```

Check filtering.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/Post_GATK_Filtering
sbatch evaluate_vcf_filtering.sh
sbatch summarize_missingness.sh
sbatch summarize_vcf_stats.sh
```

Add custom annotations for "known" vs "novel" variants and for each level of filtering ("filtered" vs "retained_filt1" vs "retained_filt2"). Then convert VCF to HDF5 format (the HDF5 format works better with scikit-allel).

- For SNPs a VCF of known variants (i.e., BOPA, 9K, 50K SNPs) was created. SNPs in the known variants VCF has the annotation "known", the remaining SNPs have the annotation "novel".
- There are too many SNPs to visualize all of them at once, so we will pull out chr1H only for visualization and do some downsampling. The purpose of visualization is to help tune filtering cutoffs and get a sense of how good of a job we are doing on filtering, so we don't need to include the full set of SNPs. One way to pull out just chr1H is using something like `bcftools view -t "chr1H_part1,chr1H_part2" raw_variants_snps.vcf > chr1H_raw_variants_snps.vcf`.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/Post_GATK_Filtering
# Prepare chr1H VCF for visualization in Jupyter Notebooks
sbatch scikit_allel_plot_prep.sh
```

Compare plots of joint annotations and unfiltered vs filtered variants using scikit-allel. This uses HDF5 (`.h5` extension) files as input. Utilizing MSI's Jupyter Notebooks since we have access to a lot more memory.

Jupyter Notebook: `Evaluate_filtering-snps.ipynb`
Run with MSI Jupyter Noteooks job profile "Mesabi High-mem - 12 cores, 128 GB, 4 hours, 180 GB local scratch", Python v3.8.3.

If phasing of some genotypes causes issues downstream, run the following to remove phase information:

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/Post_GATK_Filtering
sbatch unphase_genotypes.sh
```

---

### Methods: INDELs

First, we'll pull all indels from the raw VCF.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/Post_GATK_Filtering
sbatch select_indels.sh
# Pull out annotations and format as table for visualizing
sbatch variants_to_table_indels.sh
# Visualize annotations to get a sense of distributions for filtering
sbatch VCF_ann_viz-indels.sh
```

Filtering indels.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/Post_GATK_Filtering
./Filter_indels.sh

# Compare before and after filtering annotations
# Prepare indels VCF for visualization in Jupyter Notebooks
sbatch scikit_allel_plot_prep-indels.sh
```

Compare plots of joint annotations and unfiltered vs filtered variants using scikit-allel. This uses HDF5 (`.h5` extension) files as input. Utilizing MSI's Jupyter Notebooks since we have access to a lot more memory.

Jupyter Notebook: `Evaluate_filtering-indels.ipynb`
Run with MSI Jupyter Noteooks job profile "Mesabi High-mem - 12 cores, 128 GB, 4 hours, 180 GB local scratch", Python v3.8.3.
