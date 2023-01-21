# Analyses of 3 mutated barley lines

### VeP

Prepare vcf files if they have not been tabix indexed. VeP only works on bgzipped and tabix indexed VCF files.

Run VeP with scripts:

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/VEP
# All samples combined
./Vep_singletons-mut_snps.sh
./Vep_singletons-mut_indels.sh
# Vep for each sample individually
./Vep_singletons_by_sample-mut_snps.sh
./Vep_singletons_by_sample-mut_indels.sh

# Hybrid VCF is much larger, will need to submit as Slurm job
sbatch Vep-hybrid_snps.sh
sbatch Vep-hybrid_indels.sh
# Rare vs common SNPs
sbatch Vep-hybrid_rare_snps.sh
sbatch Vep-hybrid_common_snps.sh
# Rare vs common indels
sbatch Vep-hybrid_rare_indels.sh
sbatch Vep-hybrid_common_indels.sh

# Morex-sample2
./Vep_morex-sample2.sh
```

**Note:** When running VeP, don't use `--total_length` flag. Turning on this flag messes up the file format for BAD_Mutations `Vep_to_Subs.py` script.
