# BAD_Mutations

Run [BAD_Mutations](https://github.com/MorrellLAB/BAD_Mutations) to predict deleterious variants.

#### Step 1: Make config

```bash
# In dir: ~/Software/BAD_Mutations
./BAD_Mutations.py setup -b ~/Shared/Projects/Mutant_Barley/results/bad_mutations/cds_database -t "Hvulgare" -e 0.05 -c ~/Shared/Projects/Mutant_Barley/results/bad_mutations/config.txt
```

#### Step 2: Download CDS files

Normally, we would download the files, but since we already have the genomes available as part of another project we will use those instead. The path is: `/panfs/roc/groups/9/morrellp/shared/Projects/Selective_Sweeps/BAD_Mutations_Genomes`.

#### Step 3: Generate substitutions files

This step converts the VeP .txt files to a format that can be included in BAD_Mutations.

```bash
# In dir: ~/Software/BAD_Mutations
python ./Supporting/VeP_to_Subs.py \
       /panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/Morex_Mutants-singletons_only_missense.txt.gz \
       /panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/long_substitutions.txt \
       /panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/per-transcript_substitutions
```

#### Step 4: Generate alignments and trees

To take advantage of GNU parallel and job arrays, we will first split the *H. vulgare* all CDS fasta file into one sequence record per file.

```bash
# In dir: ~/Projects/Mutant_Barley/results/bad_mutations
~/GitHub/Barley_Mutated/02_analysis/bad_mutations/split_cds_fasta.py \
    ~/Projects/Selective_Sweeps/BAD_Mutations_Genomes/Hordeum_vulgare/Hordeum_vulgare.IBSC_v2.cds.all.fa \
    ~/Projects/Mutant_Barley/results/bad_mutations/cds_database_hvulgare
```

Now, we will generate the lists and list of lists to parallelize over.

```bash
# In dir: ~/Projects/Mutant_Barley/results/bad_mutations/cds_database_hvulgare
find $(pwd -P) -name "*.fa" | sort -V > ../align_lists/all_cds_hvulgare_list.txt
# Generate lists containing 1000 sequence records each file
split -l 1000 --numeric-suffixes all_cds_hvulgare_list.txt hvulgare_cds_list- --suffix-length=3 --additional-suffix=.txt
# Create list of lists
find $(pwd -P) -name "*list-*.txt" | sort -V > all_cds_hvulgare_list_of_lists.txt
```

Run the alignment.

```bash
sbatch --array=0-236 bad_mut_align.sh
```

