# Analyses of 3 mutated barley lines

### VeP

Prepare vcf files. VeP only works on bgzipped and tabix indexed VCF files.

```bash
# In dir: ~/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated/Filtered
module load htslib/1.9

# bgzip and index
bgzip -c mut_3_lines_filtered_singletons_only.vcf > mut_3_lines_filtered_singletons_only.vcf.gz
tabix -p vcf mut_3_lines_filtered_singletons_only.vcf.gz

bgzip -c mut_3_lines_filtered_hom_singletons_only.vcf > mut_3_lines_filtered_hom_singletons_only.vcf.gz
tabix -p vcf mut_3_lines_filtered_hom_singletons_only.vcf.gz
```

Run VeP with scripts:

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/VEP
./Vep_hom_singletons_by_sample.job
./Vep_hom_singletons_only.job
./Vep_morex-sample2.job
./Vep_singletons_by_sample.job
./Vep_singletons.job
```

**Note:** When running VeP, don't use `--total_length` flag. Turning on this flag messes up the file format for BAD_Mutations `Vep_to_Subs.py` script.
