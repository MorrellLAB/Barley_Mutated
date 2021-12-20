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
~/GitHub/Barley_Mutated/02_analysis/bad_mutations/check_align_output_file_counts.job
```

For transcripts where the `MSA_output/all_parallel_log_files/*.log` files indicate there was an error but the exit status was `0`, find the transcript name in the `MSA_output/all_parallel_log_files` list that is associated with that list number and delete the line for that transcript before re-running.

Here's an example:

The missing tree file we identified is associated with the log file `/home/morrellp/liux1299/Projects/Mutant_Barley/results/bad_mutations/MSA_output/hvulgare_cds_list-003/all_log_files/HORVU.MOREX.r2.1HG0024570.1.log`. There seems to be some error here and we want to pull it out and investigate what's going on. This is a transcript in the batch `hvulgare_cds_list-003`, so we'll go to the GNU parallel log file:

```bash
# hvulgare_cds_list-003 corresponds to index 3 in bad_mut_align.sh.3.log
vim /panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/MSA_output/all_parallel_log_files/bad_mut_align.sh.3.log
```

Then do a search for the transcript `HORVU.MOREX.r2.1HG0024570.1` within Vim and delete that line so that GNU parallel will know to re-run that transcript when we resubmit the align script. We'll repeat this process for the transcripts that did not run to completion.

For this case, we identified 11 transcripts that we want to try re-running with a higher e-value threshold. This process is a bit customized, so contact Chaochih for how she did this if it applies to your dataset.

Check to see if there are any transcripts that didn't work, if so make a note of the ones that didn't work and proceed to the predict step. It is known that some transcripts just don't work.

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
