# BAD_Mutations

Run [BAD_Mutations](https://github.com/MorrellLAB/BAD_Mutations) to predict deleterious variants.

#### Step 1: Make config

```bash
# In dir: ~/Software/BAD_Mutations
./BAD_Mutations.py setup \
    -b ~/Shared/Projects/Mutant_Barley/results/bad_mutations/cds_database \
    -t "Hvulgare" \
    -e 0.05 \
    -c ~/Shared/Projects/Mutant_Barley/results/bad_mutations/config.txt
```

#### Step 2: Download CDS files

Normally, we would download the files, but since we already have the genomes available as part of another project we will use those instead. The path is: `/panfs/roc/groups/9/morrellp/shared/Projects/Selective_Sweeps/BAD_Mutations_Genomes`.

#### Step 3: Generate substitutions files

This step converts the VeP .txt.gz files to a format that can be included in BAD_Mutations.

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations
./vep_to_subs.sh
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
# Generate lists containing 500 sequence records each file
split -l 500 --numeric-suffixes all_cds_hvulgare_list.txt hvulgare_cds_list- --suffix-length=3 --additional-suffix=.txt
# Create list of lists
find $(pwd -P) -name "*list-*.txt" | sort -V > all_cds_hvulgare_list_of_lists.txt
```

Run the alignment.

```bash
sbatch --array=0-472 bad_mut_align.sh
```

There were some transcripts that took >1 week walltime and the pasta align step still didn't finish. It is known that some transcripts just don't work, so we will make note of the 9 that didn't work and proceed to the predict step. Below is a list of ones that didn't work.

FASTA list: `/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/align_lists/hvulgare_cds_list-252.txt`

Array indices from above list: 22,27,30,32,33,34,35,36,37 corresponding to the following fasta files respectively:

```bash
/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/cds_database_hvulgare/HORVU4Hr1G053250.1.fa
/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/cds_database_hvulgare/HORVU4Hr1G053250.6.fa
/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/cds_database_hvulgare/HORVU4Hr1G053250.9.fa
/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/cds_database_hvulgare/HORVU4Hr1G053250.11.fa
/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/cds_database_hvulgare/HORVU4Hr1G053250.12.fa
/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/cds_database_hvulgare/HORVU4Hr1G053250.13.fa
/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/cds_database_hvulgare/HORVU4Hr1G053250.14.fa
/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/cds_database_hvulgare/HORVU4Hr1G053250.15.fa
/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/cds_database_hvulgare/HORVU4Hr1G053250.16.fa
```

#### Step 5: Predict substitutions

Generate a list of directories that contain the `*.subs` files. In this case, each subdirectory corresponds to one sample.

```bash
# In dir: ~/Projects/Mutant_Barley/results/bad_mutations/vep_to_subs
find $(pwd -P) -mindepth 1 -maxdepth 1 -type d > subs_dir_list.txt
```

BAD_Mutations `predict` step only works with fasta files that end in `.fasta` and NOT `.fa`. Rename fasta files.

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations
./rename_fasta_extension.sh
```

Run BAD_Mutations predict. The `.job` script stores filepaths and calls on the main script `bad_mut_predict.sh`.

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations
sbatch --array=0-2364 bad_mut_predict-M20-2-2.job
```
