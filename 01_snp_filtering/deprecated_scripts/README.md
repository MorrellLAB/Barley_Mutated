# Deprecated scripts

Scripts that were once used for a previous reference genome version but are no longer applicable for the current reference version.

---

(SKIP for Morex v3) For the phased variants vcf files relative to **Morex v1**, there were many positions with missing REF alleles (this is also a bug in 10x Genomics longranger v2.2.2). So previously, we fixed any positions with missing REF allele.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering
# IMPORTANT: this step was only necessary for when mapped to Morex v1
#qsub -t 1-4 fix_ref_allele_phased_variants.job
```

Concatenate the following VCF files output from longranger:

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering
# Mutated lines
#for i in concat_vcfs_M*.job; do qsub $i; done
# morex-sample2
#qsub concat_vcfs_morex-sample2.job
```
