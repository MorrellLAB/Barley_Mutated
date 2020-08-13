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
~/GitHub/Barley_Mutated/02_analysis/VEP/Vep_singletons.sh
~/GitHub/Barley_Mutated/02_analysis/VEP/Vep_hom_singletons_only.sh
```
