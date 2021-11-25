# BAD_Mutations

Run [BAD_Mutations](https://github.com/MorrellLAB/BAD_Mutations) to predict deleterious variants.

#### Step 1: Make config

```bash
# In dir: ~/Software/BAD_Mutations
./BAD_Mutations.py setup \
    -b ~/Shared/Projects/Mutant_Barley/results/bad_mutations/Genomes \
    -t "Hvulgare" \
    -e 0.05 \
    -c ~/Shared/Projects/Mutant_Barley/results/bad_mutations/config.txt
```

#### Step 2: Download CDS files

Download the CDS files.

```bash
~/Software/BAD_Mutations/BAD_Mutations.py fetch \
    -c /panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/config.txt

# Then remove the target species, in this case "Hvulgare"
# This is a bug that needs to be fixed, but we will manually remove it for now
cd ~/Projects/Mutant_Barley/results/bad_mutations/Genomes
rm -rf Hvulgare/
```

Make sure to exclude pangenomes and any species that gets listed twice ("duplicate" species) but are just different accessions (in this case, pick one of them to use, preferably finished assemblies). If the species gets listed twice because of female/male isolates (e.g., CpurpureusGG and CpurpureusR), for running BAD_Mutations it is better to pick the female isolate (if applicable). This information can be found through JGI Phytozome 13: https://phytozome-next.jgi.doe.gov/

```bash
# In dir: ~/Projects/Mutant_Barley/results/bad_mutations/Genomes
# Remove "duplicate" species
# We'll put it in scratch space for now in case we make a mistake
mv CpurpureusR ~/scratch/bad_mutations/duplicate_genomes
mv MpusillaCCMP ~/scratch/bad_mutations/duplicate_genomes

# Check the total number of genomes as of March 17, 2021 download
ls | wc -l
112
```

Since we had issues with runtime including all 112 species, we will reduce the number of species. But since we have already downloaded the genomes, we'll just remove the ones we don't want to include.

```bash
# In dir: ~/Projects/Mutant_Barley/results/bad_mutations
mv Genomes excluded_genomes
mkdir Genomes

# In dir: ~/Projects/Mutant_Barley/results/bad_mutations/excluded_genomes
# Only include genomes on the list
for i in $(cat ~/GitHub/Barley_Mutated/02_analysis/bad_mutations/combined_phytozome_ensembl_list.txt)
do
    # Check if directory exists
    if [[ -d "${i}" ]]; then
        mv ${i} ~/Projects/Mutant_Barley/results/bad_mutations/Genomes
    else
        echo "${i}" >> genomes_missing_from_list.txt
    fi
done
```

The following genomes on our original list wasn't already downloaded:

```bash
# In dir: ~/Projects/Mutant_Barley/results/bad_mutations/excluded_genomes
cat genomes_missing_from_list.txt 
Claxum
Platifolius
Pvaginatum
Ufusca
```

We'll exclude Claxum, Platifolius, Pvaginatum, and Ufusca since they have data restrictions.

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
module load python3/3.6.3_anaconda5.0.1

mkdir cds_database_hvulgare
~/GitHub/Barley_Mutated/02_analysis/bad_mutations/split_cds_fasta.py \
    /panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/gene_annotation/Barley_Morex_V2_gene_annotation_PGSB.HC.cds.fasta \
    ~/Projects/Mutant_Barley/results/bad_mutations/cds_database_hvulgare
```

Now, we will generate the lists and list of lists to parallelize over.

```bash
# In dir: ~/Projects/Mutant_Barley/results/bad_mutations/cds_database_hvulgare
find $(pwd -P) -name "*.fa" | sort -V > ../align_lists/all_cds_hvulgare_list.txt
# Generate lists containing 300 sequence records each file
split -l 300 --numeric-suffixes all_cds_hvulgare_list.txt hvulgare_cds_list- --suffix-length=3 --additional-suffix=.txt
# Create list of lists
find $(pwd -P) -name "*list-*.txt" | sort -V > all_cds_hvulgare_list_of_lists.txt
```

Run the alignment.

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations
sbatch --array=0-109 bad_mut_align.sh
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

To reduce the number of files we need to run, generate a list of primary transcripts only. The primary transcripts can be downloaded from Phytozome13 (https://phytozome-next.jgi.doe.gov/).

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations
./get_primary_transcripts_only.sh
```

Run BAD_Mutations predict. The `.job` script stores filepaths and calls on the main script `bad_mut_predict.sh`. The general command to submit arrays is as below, but in practice we submitted them in batches of 200 array indices for trackability.

**USEFUL:** For each batch of array indices submitted, there will be some array indices that have timed out or failed and need to be re-run. To generate a list of re-run array indices, run the script `get_re-run_array_indices.sh`.

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations
sbatch --array=0-4729 bad_mut_predict-mut_3_lines.job
```
