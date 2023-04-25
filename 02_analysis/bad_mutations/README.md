# BAD_Mutations

Run [BAD_Mutations](https://github.com/MorrellLAB/BAD_Mutations) to predict deleterious variants.

Load dependencies for BAD_Mutations.

```bash
module load python3/3.6.3_anaconda5.0.1
# Activate conda environment
source activate /home/morrellp/liux1299/.conda/envs/bad_mutations
```

## Step 1: Make config

```bash
# In dir: ~/Software/BAD_Mutations
./BAD_Mutations.py setup \
    -b ~/Shared/Projects/Mutant_Barley/results/bad_mutations/Genomes \
    -t "Hvulgare" \
    -e 0.05 \
    -c ~/Shared/Projects/Mutant_Barley/results/bad_mutations/config.txt
```

## Step 2: Download CDS files

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

## Step 3: Generate substitutions files

We may decide to go with either ANNOVAR annotations instead of VeP ones or some version of an intersect between the two.

Prepare VeP files to be converted to .subs format for BAD_Mutations. Pull out nonsynonymous variants as defined by VeP, which uses Sequence Ontology's definition (http://www.sequenceontology.org/).

This step converts the VeP .txt.gz files to a format that can be included in BAD_Mutations.

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations
./vep_to_subs-mut.sh
./vep_to_subs-hybrid_rare.sh
./vep_to_subs-hybrid_common.sh
```

## Step 4: Generate alignments and trees

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

## Step 5: Predict substitutions

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
# Full path to per transcript substitutions directory containing .subs files
#	This output is from the VeP_to_Subs.py supporting script
SUBS_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/vep_to_subs-mut/per_transcript_subs-mut_SNPs_private"
# Sample name will be used as a prefix for outputs
SAMPLE_NAME="mut_lines"
# Full path to output directory
OUT_DIR="/scratch.global/liux1299/bad_mutations/predict_output_${SAMPLE_NAME}"
# Full path to a list of primary transcripts, one per line
PRIMARY_TRANSCRIPTS="/panfs/jay/groups/9/morrellp/shared/Projects/WBDC_inversions/bad_mutations/results/phytozome13_download_V3_primary_transcript/hvulgare_primary_transcripts_only.txt"

# Prepare subs list that intersects with primary transcripts list
./intersect_primary_transcripts_and_subs.sh ${SUBS_DIR} ${OUT_DIR} ${PRIMARY_TRANSCRIPTS}

# Run bad mutations predict
sbatch --array=0-209 bad_mut_predict-mut.sh
```

Check which subdirectories have problematic predictions (i.e., error when this step was run resulting in output files with error messages).

```bash
# In dir: ~/scratch/bad_mutations/predict_output_mut_lines
for i in $(find hvulgare_cds_list-* -mindepth 1 -maxdepth 1 -not -empty -type d); do
    ls $i/*
done
# Output
# hvulgare_cds_list-015/problematic_predictions/HORVU.MOREX.r3.1HG0058740.1_Predictions.txt
# hvulgare_cds_list-078/problematic_predictions/HORVU.MOREX.r3.3HG0305520.1_Predictions.txt
# hvulgare_cds_list-110/problematic_predictions/HORVU.MOREX.r3.5HG0428210.1_Predictions.txt
# hvulgare_cds_list-145/problematic_predictions/HORVU.MOREX.r3.6HG0565340.1_Predictions.txt
# hvulgare_cds_list-179/problematic_predictions/HORVU.MOREX.r3.7HG0696890.1_Predictions.txt

# Error types
# 1) Internal error, dumping the offending likelihood function to /tmp/hyphy.dump 
# 2) The leaf of the tree:givenTree labeled _41562 had no match in the data. Please make sure that all leaf names correspond to a sequence name in the data file.
```

Repeat for Hybrid rare and Hybrid common SNPs.

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations
# Full path to per transcript substitutions directory containing .subs files
#	This output is from the VeP_to_Subs.py supporting script
SUBS_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/vep_to_subs-hybrid13/per_transcript_subs-hybrid_SNPs_rare"
# Sample name will be used as a prefix for outputs
SAMPLE_NAME="hybrid_rare"
# Full path to output directory
OUT_DIR="/scratch.global/liux1299/bad_mutations/predict_output_${SAMPLE_NAME}"
# Full path to a list of primary transcripts, one per line
PRIMARY_TRANSCRIPTS="/panfs/jay/groups/9/morrellp/shared/Projects/WBDC_inversions/bad_mutations/results/phytozome13_download_V3_primary_transcript/hvulgare_primary_transcripts_only.txt"

# Prepare subs list that intersects with primary transcripts list
./intersect_primary_transcripts_and_subs.sh ${SUBS_DIR} ${OUT_DIR} ${PRIMARY_TRANSCRIPTS}

# Run bad mutations predict
sbatch --array=0-209 bad_mut_predict-hybrid_rare.sh
# Continue running timeout indices
sbatch --array=0-4,8,10,14,16-17,19-35,37-38,42-111,113-193 bad_mut_predict-hybrid_rare.sh
# Continue running timeout indices while troubleshooting failed indices
sbatch --array=0-4,16,20-31,45-48,50-62,64,69,72,75-83,85-88,90-91,98,100-102,105-107,109-110,119,124,126-128,132,135-140,159-163,165,167-169,173,187-188,191 bad_mut_predict-hybrid_rare.sh
# Re-run failed indices after update to code
sbatch --array=17,19,32-35,37-38,42-44,49,63,65-68,70-71,73-74,84,89,92-97,99,103-104,108,111,113-118,120-123,125,129-131,133-134,141-158,164,166,170-172,174-186,189-190,192-193 bad_mut_predict-hybrid_rare.sh
```

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations
# Full path to per transcript substitutions directory containing .subs files
#	This output is from the VeP_to_Subs.py supporting script
SUBS_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/vep_to_subs-hybrid13/per_transcript_subs-hybrid_SNPs_common"
# Sample name will be used as a prefix for outputs
SAMPLE_NAME="hybrid_common"
# Full path to output directory
OUT_DIR="/scratch.global/liux1299/bad_mutations/predict_output_${SAMPLE_NAME}"
# Full path to a list of primary transcripts, one per line
PRIMARY_TRANSCRIPTS="/panfs/jay/groups/9/morrellp/shared/Projects/WBDC_inversions/bad_mutations/results/phytozome13_download_V3_primary_transcript/hvulgare_primary_transcripts_only.txt"

# Prepare subs list that intersects with primary transcripts list
./intersect_primary_transcripts_and_subs.sh ${SUBS_DIR} ${OUT_DIR} ${PRIMARY_TRANSCRIPTS}

# Run bad mutations predict
sbatch --array=0-209 bad_mut_predict-hybrid_common.sh
# Continue running timeout indices
sbatch --array=0-5,11-17,19-32,34-36,41-61,64-65,71-90,99,101-111,113-114,118-119,122-142,147-148,152-154,156-174,176,181-184,186-193 bad_mut_predict-hybrid_common.sh
# Re-run failed and debug
sbatch --array=100 bad_mut_predict-hybrid_common.sh
# Continue running timeout indices
sbatch --array=0,24,27,35,43-48,50-61,76-78,80-83,85,88,102-103,106-107,109,111,123,125-132,135-142,157,159-165,167-169,171,184,190-191,193 bad_mut_predict-hybrid_common.sh
# Continue running timeout indices
sbatch --array=142 bad_mut_predict-hybrid_common.sh
```

## Step 6: Compile predictions

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations
sbatch bad_mut_compile_predict-mut.sh
sbatch bad_mut_compile_predict-hyb_rare.sh
sbatch bad_mut_compile_predict-hyb_common.sh
```

## Step 7: Post processing and Visualization

**Get Deleterious vs Tolerated:** Merge BAD_Mutations predict compiled report and VeP report. Then, annotate as "Deleterious" vs "Tolerated".

Mutated lines

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations

# Mutated
# Criteria for annotating as "Deleterious" vs "Tolerated"
MIN_SEQ="10"
MAX_CONSTRAINT="1"
P_CUTOFF="0.05"
# Note: this should be the number of nonsynonymous SNPs given as input to the BAD_Mutations predict step
#   and NOT the number of lines in the compiled report
# Can get number from the primary transcript names intersect input file to BAD_Mutations predict
# /scratch.global/liux1299/bad_mutations/predict_output_mut_lines/primary_transcript_intersect_names_only.txt
NUM_CODONS_TESTED_MUT="611"

./prep_del_vs_tol_vep.sh \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/predictions/mut_snps_private_Combined_Report.txt \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/HC_LC_gff_SNPs_private_all_samples/mut8_and_3mut10xGenomics.SNPs.private.txt \
    mut_snps_private \
    ${MIN_SEQ} ${MAX_CONSTRAINT} ${P_CUTOFF} ${NUM_CODONS_TESTED_MUT}
```

Repeat above for hybrid rare and common nonsynonymous SNPs.

Rare

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations

# Rare
# Criteria for annotating as "Deleterious" vs "Tolerated"
MIN_SEQ="10"
MAX_CONSTRAINT="1"
P_CUTOFF="0.05"
# Note: this should be the number of nonsynonymous SNPs given as input to the BAD_Mutations predict step
#   and NOT the number of lines in the compiled report
# Can get number from the primary transcript names intersect input file to BAD_Mutations predict
# /scratch.global/liux1299/bad_mutations/predict_output_hybrid_rare/primary_transcript_intersect_names_only.txt
NUM_CODONS_TESTED_RARE="9716"

./prep_del_vs_tol_vep.sh \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/predictions/hybrid_rare_Combined_Report.txt \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/HC_LC_gff_SNPs_hybrid13_rare/hybrid13.SNPs.rare.txt \
    hybrid13_rare_snps \
    ${MIN_SEQ} ${MAX_CONSTRAINT} ${P_CUTOFF} ${NUM_CODONS_TESTED_RARE}
```

Common

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations

# Common
# Criteria for annotating as "Deleterious" vs "Tolerated"
MIN_SEQ="10"
MAX_CONSTRAINT="1"
P_CUTOFF="0.05"
# Note: this should be the number of nonsynonymous SNPs given as input to the BAD_Mutations predict step
#   and NOT the number of lines in the compiled report
# Can get number from the primary transcript names intersect input file to BAD_Mutations predict
# /scratch.global/liux1299/bad_mutations/predict_output_hybrid_common/primary_transcript_intersect_names_only.txt
NUM_CODONS_TESTED_COMMON="14537"

./prep_del_vs_tol_vep.sh \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/predictions/hybrid_common_Combined_Report.txt \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/HC_LC_gff_SNPs_hybrid13_common/hybrid13.SNPs.common.txt \
    hybrid13_common_snps \
    ${MIN_SEQ} ${MAX_CONSTRAINT} ${P_CUTOFF} ${NUM_CODONS_TESTED_COMMON}
```

#### Prepare files for plotting dSNPs per codon

Count number of codons per 10 Mbp windows.

```bash
GENOME_FILE="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta.fai"
# Positions of transcripts from GFF
GFF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/annotation/HvulgareMorex_702_V3.gene.gff3"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/visualization"
WIN_SIZE="10000000"

~/GitHub/Barley_Mutated/02_analysis/bad_mutations/num_codons_per_window.sh ${GENOME_FILE} ${GFF} ${OUT_DIR} ${WIN_SIZE} "CDS"

~/GitHub/Barley_Mutated/02_analysis/bad_mutations/num_codons_per_window.sh ${GENOME_FILE} ${GFF} ${OUT_DIR} ${WIN_SIZE} "mRNA"
```

Mutated lines dataset.

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations
./prep_dsnps_per_codon.sh \
    "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/predictions/mut_snps_private_deleterious_vs_tolerated.txt" \
    "/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta.fai" \
    "/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/annotation/HvulgareMorex_702_V3.gene.gff3" \
    "10000000" \
    "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/visualization" \
    "mut_snps_private"
```

Rare:

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations
./prep_dsnps_per_codon.sh \
    "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/predictions/hybrid13_rare_snps_deleterious_vs_tolerated.txt" \
    "/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta.fai" \
    "/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/annotation/HvulgareMorex_702_V3.gene.gff3" \
    "10000000" \
    "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/visualization" \
    "hybrid13_rare_snps"
```

Common:

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations
./prep_dsnps_per_codon.sh \
    "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/predictions/hybrid13_common_snps_deleterious_vs_tolerated.txt" \
    "/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta.fai" \
    "/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/annotation/HvulgareMorex_702_V3.gene.gff3" \
    "10000000" \
    "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/visualization" \
    "hybrid13_common_snps"
```

#### Count SNPs: Deleterious vs Tolerated and Nonsynonymous vs Synonymous

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations
module load python3/3.8.3_anaconda2020.07_mamba

# Mutated
./count_dSNP_per_individual.py \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.SNPs.private.vcf.gz \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/HC_LC_gff_SNPs_private_all_samples/mut8_and_3mut10xGenomics.SNPs.private.txt \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/predictions/mut_snps_private_deleterious_vs_tolerated.txt \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/visualization \
    mut_private

# Rare
# Run as interactive job
srun -N 1 --ntasks-per-node=8 --mem=36gb --tmp=22gb -t 3:00:00 -p interactive --pty bash
module load python3/3.8.3_anaconda2020.07_mamba
./count_dSNP_per_individual.py \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/hybrid13.SNPs.rare.vcf.gz \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/HC_LC_gff_SNPs_hybrid13_rare/hybrid13.SNPs.rare.txt \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/predictions/hybrid13_rare_snps_deleterious_vs_tolerated.txt \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/visualization \
    hybrid13.SNPs.rare

# Common
# Run as interactive job
srun -N 1 --ntasks-per-node=8 --mem=36gb --tmp=22gb -t 3:00:00 -p interactive --pty bash
module load python3/3.8.3_anaconda2020.07_mamba
./count_dSNP_per_individual.py \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/hybrid13.SNPs.common.vcf.gz \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/HC_LC_gff_SNPs_hybrid13_common/hybrid13.SNPs.common.txt \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/predictions/hybrid13_common_snps_deleterious_vs_tolerated.txt \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/visualization \
    hybrid13.SNPs.common
```

Count indels (even though we didn't run BAD_Mutations on indels since it takes SNPs, still wanted to get per sample counts with the same counting approach as in the count_dSNP_per_individual.py script).

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/bad_mutations
module load python3/3.8.3_anaconda2020.07_mamba

# Mutated
./count_indels_per_individual.py \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.INDELs.private.vcf.gz \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/visualization \
    mut_private_indels

# Rare
./count_indels_per_individual.py \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/hybrid13.INDELs.rare.vcf.gz \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/visualization \
    hybrid_rare_indels

# Common
./count_indels_per_individual.py \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/hybrid_rare_vcfs/hybrid13.INDELs.common.vcf.gz \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/visualization \
    hybrid_common_indels
```

See `03_figures` directory for plottings scripts.

## Unused steps during early exploration

Prepare VeP reports for generating figures of BAD_Mutations results.

```bash
# Dependencies
module load perl/modules.centos7.5.26.1
module load htslib/1.9
export PATH=$PATH:/panfs/jay/groups/9/morrellp/shared/Software/ensembl-vep-release-108.1
```

Pull out synonymous variants.

```bash
MUT_VEP_REPORT="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/HC_LC_gff_SNPs_private_all_samples/mut8_and_3mut10xGenomics.SNPs.private.txt"

MUT_PREFIX=$(basename $MUT_VEP_REPORT .txt)
MUT_OUT_DIR=$(dirname $MUT_VEP_REPORT)

filter_vep -i $MUT_VEP_REPORT -o $MUT_OUT_DIR/$MUT_PREFIX.synonymous.txt -filter "Consequence is synonymous_variant"
# Gzip output file
gzip $MUT_OUT_DIR/$MUT_PREFIX.synonymous.txt
```

Get the total number of codons for SNPs per codon figure.

```bash
module load samtools/1.9

# Use primary transcripts in cds from JGI
fasta_file="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/annotation/HvulgareMorex_702_V3.cds_primaryTranscriptOnly.fa"

samtools faidx $fasta_file > /panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/annotation/HvulgareMorex_702_V3.cds_primaryTranscriptOnly.fa.fai
```

The length of each cds primary transcript is in the file:

```bash
CDS_PT_FAI="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/annotation/HvulgareMorex_702_V3.cds_primaryTranscriptOnly.fa.fai"
# This list was used to determine primary transcripts to run for BAD_Mutations
# Same transcripts as in the CDS_PT_FAI file
CDS_PT_NAMES="/panfs/jay/groups/9/morrellp/shared/Projects/WBDC_inversions/bad_mutations/results/phytozome13_download_V3_primary_transcript/hvulgare_primary_transcripts_only.txt"
```

Generate intervals for pseudomolecules positions from `.fai` file.

```bash
module load bedtools/2.29.2

GENOME_FILE="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta.fai"
# Window size in bp
WIN_SIZE="10000000"
TMP_DIR=~/scratch/temp_visualization

# In dir: ~/Projects/Mutant_Barley/results/bad_mutations/visualization
# Create non-overlapping windows
bedtools makewindows -g ${GENOME_FILE} -w ${WIN_SIZE} | grep -v "chrUn" > Barley_MorexV3_pseudomolecules.10Mb_windows.bed
# Split file into one line per file
split -l 1 --numeric-suffixes --additional-suffix=".bed" --suffix-length=3 Barley_MorexV3_pseudomolecules.10Mb_windows.bed ${TMP_DIR}/Barley_MorexV3_pseudomolecules.10Mb_windows.
```
