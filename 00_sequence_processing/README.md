# Sequence Processing

Scripts used for sequence processing are stored in subdirectories and categorized by sequencing technology.

- `10x_Genomics_Morex_v3`
- `Nanopore`
- `PacBio`

Config files used with [`sequence_handling`](https://github.com/MorrellLAB/sequence_handling) for processing Illumina WGS sequencing are stored in this directory.

To generate the chromosome intervals to parallelize over, I ran the following:

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing
module load python3/3.7.4_anaconda2019.10
# Split chromosomes into smaller regions for GenomicsDBImport parallelization
./split_chrom_fai.py /panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta.fai 2 > ~/Alignments/mut8_and_hybrid_barley/barley_chromosomes_split.intervals
```
