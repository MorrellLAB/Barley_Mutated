# BAD_Mutations

Run [BAD_Mutations](https://github.com/MorrellLAB/BAD_Mutations) to predict deleterious variants.

Load dependencies for BAD_Mutations.

```bash
module load python3/3.6.3_anaconda5.0.1
# Activate conda environment
source activate /home/morrellp/liux1299/.conda/envs/bad_mutations
```

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

*Note:* If you have a list of favorite species you have picked out, you can modify the `BAD_Mutations/lrt_predict/Fetch/phytozome_species.py` script to only include your favorite species. Then, only your favorite species will be downloaded. We'll use various shortcuts in VSCode (including adding text to the start and end of every line: https://shouts.dev/how-to-add-text-at-start-and-end-of-all-lines-in-vs-code) to format our `combined_phytozome_ensembl_list.txt` to match the syntax in `BAD_Mutations/lrt_predict/Fetch/phytozome_species.py`. Remember, this is the phytozome list only not the ensemble list, so we'll need to remove the lines that contain species naming from Ensembl. These should just be `Aegilops_tauschii`, `Leersia_perrieri`, and `Triticum_urartu` on our combined list.

```bash
~/Software/BAD_Mutations/BAD_Mutations.py fetch \
    -c /panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/config.txt

# Then remove the target species, in this case "Hvulgare"
# This is a bug that needs to be fixed, but we will manually remove it for now
cd ~/Projects/Mutant_Barley/results/bad_mutations/Genomes
rm -rf Hvulgare/
```

Since we had issues with runtime when including all 112 species available for download at the time, we will reduce the number of species. But since we have already downloaded the genomes, we'll just keep the ones we want to include and remove the ones we don't want to include.

We selected only Angiosperm genomes available on Phytozome 13 (https://phytozome-next.jgi.doe.gov/) because more closely related species adds to the quality of annotation more than very distantly related species. Genomes that aren't Angiosperms also seem to be very difficult to align, which we suspect was causing the extremely long runtimes (some genes took longer than 10 days to align).

```bash
# In dir: ~/Projects/Mutant_Barley/results/bad_mutations
mkdir excluded_genomes
#mv Genomes/* excluded_genomes/

# In dir: ~/Projects/Mutant_Barley/results/bad_mutations/Genomes
cd Genomes
# Only include genomes on the list
for i in $(ls -d * | grep -v ".txt")
do
    # If genome is not on our list, move it to the excluded_genomes directory
    if ! grep -w -q "${i}" ~/GitHub/Barley_Mutated/02_analysis/bad_mutations/combined_phytozome_ensembl_list.txt
    then
        # Move excluded genome
        echo "Moving excluded genome: ${i}"
        mv "${i}" ~/Projects/Mutant_Barley/results/bad_mutations/excluded_genomes
    fi
done

# Check which genomes on our list was not downloaded by BAD_Mutations fetch
for i in $(cat ~/GitHub/Barley_Mutated/02_analysis/bad_mutations/combined_phytozome_ensembl_list.txt)
do
    if ! [[ -d "${i}" ]]; then
        echo "${i}" >> genomes_not_downloaded.txt
    fi
done
```

The following genomes on our original list wasn't already downloaded:

```bash
# In dir: ~/Projects/Mutant_Barley/results/bad_mutations/excluded_genomes
cat genomes_not_downloaded.txt 
Claxum
Pacutifolius
Platifolius
Pvaginatum
Ufusca
```

We'll add these 5 species manually to the BAD_Mutations script that contains a list of Phytozome species. IMPORTANT!!! Make sure to remove species not on our list too! Otherwise, the species not on our list will also be downloaded too. We'll use various shortcuts in VSCode (including adding text to the start and end of every line: https://shouts.dev/how-to-add-text-at-start-and-end-of-all-lines-in-vs-code) to format our `combined_phytozome_ensembl_list.txt` to match the syntax in `BAD_Mutations/lrt_predict/Fetch/phytozome_species.py`. Remember, this is the phytozome list only not the ensemble list, so we'll need to remove the lines that contain species naming from Ensembl. These should just be `Aegilops_tauschii`, `Leersia_perrieri`, and `Triticum_urartu` on our combined list.

```bash
# In MSI dir: ~/Software/BAD_Mutations/lrt_predict/Fetch
# Create a copy of the phytozome list that was originally in the directory in case something goes wrong
cp phytozome_species.py original_copy_phytozome_species.py
# Copy our custom list of species and name it the same way as the original file
cp ~/GitHub/Barley_Mutated/02_analysis/bad_mutations/my_phytozome_species.py phytozome_species.py

# Note: Now looking back, I could have just generated a version of the phytozome_species.py
#   that includes just my favorite species, but I already had all 112 species downloaded so
#   I used the code above to remove species we wanted to exclude.

# Now, re-run the fetch command to fetch remaining genomes
~/Software/BAD_Mutations/BAD_Mutations.py fetch \
    -c /panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/config.txt

# One last check to see if we are still missing some genomes we wanted to include
for i in $(cat ~/GitHub/Barley_Mutated/02_analysis/bad_mutations/combined_phytozome_ensembl_list.txt)
do
    if ! [[ -d "${i}" ]]; then
        echo "${i}"
    fi
done
```

We'll also need to check that each species directory doesn't contain more than two versions of the genome and only contains the latest version of the genome.

We'll proceed with 72 species for alignment located in the following directory:

```bash
/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/Genomes
```

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
