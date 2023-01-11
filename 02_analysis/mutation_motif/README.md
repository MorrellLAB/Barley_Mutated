# Context of Mutations analysis

Run [Mutation Motif](https://github.com/HuttleyLab/MutationMotif/) to analyze the context of mutations.

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/mutation_motif
# Prepare counts table format
sbatch prep_counts_table-mut_snps.sh

# Prepare list of counts separated by mutation direction
# There should be 12 files total
cd ~/Projects/Mutant_Barley/results/mutation_motif/mut8_and_3mut10xGenomics.SNPs.private/counts_tables
realpath *to*.txt | sort -uV > separate_counts_file_list.txt

# In dir: ~/GitHub/Barley_Mutated/02_analysis/mutation_motif
# Run analyses and generate plots
./context_analysis-mut_snps.sh
```

Repeat for hybrid rare SNPs.

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/mutation_motif
# Prepare counts table format
sbatch prep_counts_table-hybrid_rare_snps.sh

# Prepare list of counts separated by mutation direction
# There should be 12 files total
cd ~/Projects/Mutant_Barley/results/mutation_motif/hybrid13.SNPs.rare/counts_tables
realpath *to*.txt | sort -uV > separate_counts_file_list.txt

# In dir: ~/GitHub/Barley_Mutated/02_analysis/mutation_motif
# Run analyses and generate plots
./context_analysis-hybrid_rare_snps.sh
```

Repeat for hybrid common SNPs.

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/mutation_motif
# Prepare counts table format
sbatch prep_counts_table-hybrid_common_snps.sh

# Prepare list of counts separated by mutation direction
# There should be 12 files total
cd ~/Projects/Mutant_Barley/results/mutation_motif/hybrid13.SNPs.common/counts_tables
realpath *to*.txt | sort -uV > separate_counts_file_list.txt

# Run analyses and generate plots
./context_analysis-hybrid_common_snps.sh
```
