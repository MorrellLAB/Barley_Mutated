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

We may decide to go with either ANNOVAR annotations instead of VeP ones or some version of an intersect between the two.

Prepare VeP files to be converted to .subs format for BAD_Mutations. Pull out nonsynonymous variants as defined by VeP, which uses Sequence Ontology's definition (http://www.sequenceontology.org/).

This step converts the VeP .txt.gz files to a format that can be included in BAD_Mutations.

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations
./vep_to_subs.sh
```

#### Step 4: Generate alignments and trees

> **IMPORTANT NOTE:** This is how we originally ran BAD_Mutations relative to Morex_v2 but then we decided to move to Morex v3. For Morex_v3, Giulia ran this step and generated the alignments and trees. We utilized Giulia's output from this step to run predictions.

To take advantage of GNU parallel and job arrays, we will first split the *H. vulgare* all CDS fasta file into one sequence record per file.

```bash
# In dir: ~/Projects/Mutant_Barley/results/bad_mutations
module load python3/3.6.3_anaconda5.0.1

mkdir cds_database_hvulgare
~/GitHub/Barley_Mutated/02_analysis/bad_mutations/split_cds_fasta.py \
    /panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/gene_annotation/Barley_Morex_V2_gene_annotation_PGSB.all.cds.fasta \
    ~/Projects/Mutant_Barley/results/bad_mutations/cds_database_hvulgare
```

Now, we will generate the lists and list of lists to parallelize over.

```bash
# In dir: ~/Projects/Mutant_Barley/results/bad_mutations/cds_database_hvulgare
find $(pwd -P) -name "*.fa" | sort -V > ../align_lists/all_cds_hvulgare_list.txt

# In dir: ~/Shared/Projects/Mutant_Barley/results/bad_mutations/align_lists
# Generate lists containing 400 sequence records each file
split -l 400 --numeric-suffixes all_cds_hvulgare_list.txt hvulgare_cds_list- --suffix-length=3 --additional-suffix=.txt
# Create list of lists
find $(pwd -P) -name "*list-*.txt" | sort -V > all_cds_hvulgare_list_of_lists.txt
```

Run the alignment. Our max array index should be the last list number is our lists of lists above (this is because our list numbering starts at "000"), but the max array index will also get printed to the `*.out` file that gets generated once the job starts running.

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations
# For barley Morex v2, a set of 400 transcripts in each job array took less than 8 hours
sbatch --array=0-159 bad_mut_align.sh
```

After the first pass of running all array indices is done, check if there are any array indices that timed out.

```bash
# The "9137097" is the job ID of the Slurm job submitted above
~/GitHub/Barley_Mutated/02_analysis/bad_mutations/get_re-run_array_indices.sh 9137097
# Output
4,11,22-23,37,49,52,57,60,71,76,79,82-83,87,89,94,105,114,120,125
```

The re-run array indices above are formatted correctly for re-submitting some array indices.

```bash
sbatch --array=4,11,22-23,37,49,52,57,60,71,76,79,82-83,87,89,94,105,114,120,125 bad_mut_align.sh
```

There were still some that timed out. Repeat the above.

```bash
~/GitHub/Barley_Mutated/02_analysis/bad_mutations/get_re-run_array_indices.sh 9218771
4,22-23,37,60,71,79,89,105,114

sbatch --array=4,22-23,37,60,71,79,89,105,114 bad_mut_align.sh
```

There were still some that timed out, based on the output messages, it looks like they just need longer walltime. We'll increase to 48 hours.

```bash
~/GitHub/Barley_Mutated/02_analysis/bad_mutations/get_re-run_array_indices.sh 9255778
4,23,37,60,71,105,114

sbatch --array=4,23,37,60,71,105,114 bad_mut_align.sh
```

Again, there were five genes that timed out. We'll increase walltime to 90 hours.

```bash
~/GitHub/Barley_Mutated/02_analysis/bad_mutations/get_re-run_array_indices.sh 9292552
4,37,71,105,114

sbatch --array=4,37,71,105,114 bad_mut_align.sh
```

A few are taking longer than 90 hours.

```bash
~/GitHub/Barley_Mutated/02_analysis/bad_mutations/get_re-run_array_indices.sh 9425464
4,37,105

# In dir: ~/Projects/Mutant_Barley/results/bad_mutations/MSA_output
# Check which transcripts these are
cat temp_msa_output_other_error_log_files.txt
# Output
/home/morrellp/liux1299/Projects/Mutant_Barley/results/bad_mutations/MSA_output/hvulgare_cds_list-004/all_log_files/HORVU.MOREX.r2.1HG0018220.1.log
/home/morrellp/liux1299/Projects/Mutant_Barley/results/bad_mutations/MSA_output/hvulgare_cds_list-037/all_log_files/HORVU.MOREX.r2.2HG0150370.1.log
/home/morrellp/liux1299/Projects/Mutant_Barley/results/bad_mutations/MSA_output/hvulgare_cds_list-105/all_log_files/HORVU.MOREX.r2.5HG0421780.1.log

# In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations
# Based on the log files, these just need longer to run
# We'll try increasing the walltime to 240 hours just for these 3
sbatch --array=4,37,105 bad_mut_align.sh
```

Two are taking longer than 240 hours.

```bash
~/GitHub/Barley_Mutated/02_analysis/bad_mutations/get_re-run_array_indices.sh 9660527
4,105

# In dir: ~/Projects/Mutant_Barley/results/bad_mutations/MSA_output
# Check which transcripts these are
cat temp_msa_output_other_error_log_files.txt
# Output
/home/morrellp/liux1299/Projects/Mutant_Barley/results/bad_mutations/MSA_output/hvulgare_cds_list-004/all_log_files/HORVU.MOREX.r2.1HG0018220.1.log
/home/morrellp/liux1299/Projects/Mutant_Barley/results/bad_mutations/MSA_output/hvulgare_cds_list-105/all_log_files/HORVU.MOREX.r2.5HG0421780.1.log

# In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations
# Based on the log files, they just need longer to run
# We'll try increasing the walltime to 480 hours for the last 2
sbatch --array=4,105 bad_mut_align.sh
```

*Note:* BAD_Mutations align combines both Slurm job arrays and GNU parallel. This allows re-submitting the same array index and picking up where the job left off if we run out of walltime. It does this by keeping a GNU parallel log file in the `${OUT_DIR}/all_parallel_log_files`. Each array index will have its own log file that track the exit status of each GNU parallel task.

Each subdirectory in `MSA_Output` a subdirectory called `all_log_files` (e.g., `MSA_Output/hvulgare_cds_list-000/all_log_files`, `MSA_Output/hvulgare_cds_list-001/all_log_files`, `MSA_Output/hvulgare_cds_list-002/all_log_files`, etc.) that keeps a log of the align output (i.e., the stdout from BAD_Mutations align) for each transcript in that batch. The log files have basenames that match the name of the transcript in that batch. For example `HORVU.MOREX.r2.1HG0000020.1.log` is the log file for a transcript from the list `/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/align_lists/hvulgare_cds_list-000.txt`. This makes it a little easier to troubleshoot if needed. The `MSA_output/all_parallel_log_files` is a log to keep track of the parallel tasks being run and where to pick up the run in case we run out of walltime and need to re-submit the job.

Check that we have the expected number of `*.fa` and `*.tree` files written to the `MSA_Output` directory. These files occur in pairs, so if each batch has 400 transcripts run, we would expect 400 `*.fa` files and 400 `*.tree` files in the output directory. CAUTION: Exit statuses shouldn't be the sole metric you use to determine if a job ran to completion successfully, which is why we are running these additional checks.

```bash
# In dir: ~/Projects/Mutant_Barley/results/bad_mutations/MSA_output
# First pass quick check
for i in $(ls -d hvulgare_cds_list-*)
do
    echo $i
    ls ${i}/*.tree | wc -l
done

# The number of *.tree files in each subdirectory should be the same as the
#   input list of lists
wc -l ~/Projects/Mutant_Barley/results/bad_mutations/align_lists/hvulgare_cds_list-*.txt

# We'll programmatically check that we have the expected number of output files
~/GitHub/Barley_Mutated/02_analysis/bad_mutations/check_align_output_file_counts.sh \
    /panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/align_lists/all_cds_hvulgare_list_of_lists.txt \
    ~/Projects/Mutant_Barley/results/bad_mutations/MSA_output
```

This check outputs of the following files in the `MSA_output` directory:

```bash
temp_msa_output_evalue_error_log_files.txt
temp_msa_output_problem_dirs.txt
temp_msa_output_missing_transcripts_log_files.txt
temp_msa_output_unexpected_symbol_error_log_files.txt
temp_msa_output_other_error_log_files.txt
```

The `temp_msa_output_evalue_error_log_files.txt` contains a list of transcript log files with the following error:

```
CRITICAL	Could not find any BLAST hits! Try raising the E-value threshold for homology.
```

There were 203 CDS sequences with this type of e-value threshold error. We can exclude these from the predictions because they are likely annotation errors in the barley genome or are genes that couldn't be anntoated in related species and are missing from other plant peptide sequences.

Check to see if there are any other transcripts that didn't work, if so make a note of the ones that didn't work and proceed to the predict step. It is known that some transcripts just don't work.

We'll need to regenerate the `FASTA_LIST_OF_LISTS` to exclude problematic CDS sequences so our predict script works as intended.didn't work and proceed to the predict step. It is known that some transcripts just don't work.

#### Step 5: Predict substitutions

The following parts inside the block quote was already run by Giulia since we used her align and trees output:

> Generate a list of directories that contain the `*.subs` files. In this case, each subdirectory corresponds to one sample.

> ```bash
> # In dir: ~/Projects/Mutant_Barley/results/bad_mutations/vep_to_subs
> find $(pwd -P) -mindepth 1 -maxdepth 1 -type d > subs_dir_list.txt
> ```

> BAD_Mutations `predict` step only works with fasta files that end in `.fasta` and NOT `.fa`. Rename fasta files.

> ```bash
> # In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations
> ./rename_fasta_extension.sh
> ```

> To reduce the number of files we need to run, generate a list of primary transcripts only. The primary transcripts can be downloaded from Phytozome13 (https://phytozome-next.jgi.doe.gov/).

> ```bash
> # In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations
> ./get_primary_transcripts_only.sh
> ```

Run BAD_Mutations predict. The `bad_mut_predict-mut.sh` script stores filepaths and calls on the main script `bad_mut_predict.sh`. The general command to submit arrays is as below, but in practice we submitted them in batches of 200 array indices for trackability.

**USEFUL:** For each batch of array indices submitted, there will be some array indices that have timed out or failed and need to be re-run. To generate a list of re-run array indices, run the script `get_re-run_array_indices.sh`.

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations
sbatch --array=0-209 bad_mut_predict-mut.sh
```

#### Step 6: Compile predictions


