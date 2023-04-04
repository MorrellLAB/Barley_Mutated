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

---

### Motif frequency

Prepare fasta with uncallable regions masked.

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/mutation_motif
./mask_fasta_uncallable.sh
```

Use [EMOSS compseq](https://www.bioinformatics.nl/cgi-bin/emboss/help/compseq) to calculate motif frequency.

```bash
# In dir: ~/Projects/Mutant_Barley/results/emboss_compseq
module load emboss/6.6.0
# dinucleotides
compseq -sequence Barley_MorexV3_pseudomolecules_parts.masked_uncallable.fasta -word 2 -outfile dinucleotide_morex_v3.masked_uncallable.comp
# trinucleotides
compseq -sequence Barley_MorexV3_pseudomolecules_parts.masked_uncallable.fasta -word 3 -outfile trinucleotide_morex_v3.masked_uncallable.comp
```