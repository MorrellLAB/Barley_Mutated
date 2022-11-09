# Context of Mutations analysis

Run [Mutation Motif](https://github.com/HuttleyLab/MutationMotif/) to analyze the context of mutations.

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/mutation_motif
# Prepare counts table format
sbatch prep_counts_table-mut_snps.sh
# Run analyses and generate plots
./context_analysis-mut_snps.sh
```

Repeat for hybrid rare SNPs.

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/mutation_motif
# Prepare counts table format
sbatch prep_counts_table-hybrid_rare_snps.sh
# Run analyses and generate plots
./context_analysis-hybrid_rare_snps.sh
```

Repeat for hybrid common SNPs.

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/mutation_motif
# Prepare counts table format
sbatch prep_counts_table-hybrid_common_snps.sh
# Run analyses and generate plots
./context_analysis-hybrid_common_snps.sh
```
