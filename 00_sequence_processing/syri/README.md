# Minimap2 alignment and Syri

Get a sense of the amount of variation present between replicates used for the Morex reference genome assembly process. This should give us a rough baseline of how many differences to expect between the morex sample sequenced with 10x Genomics in our study and the morex reference.

Prepare FASTA files split by chromosome for faster processing and easier visualization/comparison downstream using syri.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/syri
sbatch split_fasta_by_chr.sh
```

Run minimap2 for alignment.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/syri
sbatch --array=0-6 minimap2_asm5_morex.sh
```

Run Syri using PAF files as input.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/syri
sbatch --array=0-6 run_syri.sh
```

Plot syri output using [plotsr](https://github.com/schneebergerlab/plotsr).

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/syri
sbatch run_plotsr.sh
```
