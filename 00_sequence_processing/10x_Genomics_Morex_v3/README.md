# WGS 10x Genomics Data Mutated Lines

## Longranger

[`longranger wgs` documentation](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/using/wgs)

#### Step 1: Prepare reference fasta file

Modify and run `make_ref.sh` script using `sbatch` command to submit the script as a Slurm job to MSI.

#### Step 2: Run Longranger wgs mode

Modify the user email and the full path to the directory containing 10x Genomics compatible reference directory in the following three scripts. The remaining input arguments can be left as is.

- `run_longranger_wgs-M01.sh`
- `run_longranger_wgs-M20.sh`
- `run_longranger_wgs-M29.sh`
- `run_longranger_wgs-Morex.sh`

After modifying the scripts, use `sbatch` to submit each script. If you exceed the walltime after the first submission, re-submit the exact same script again and keep repeating this until the analysis finishes. 10x Genomics `longranger` has a feature where it can pick up processing from an "orphaned state". Previously, it has taken roughly 1 month to complete a run for barley samples.

#### Step 3: Fix output filenames

First, we'll add the sample name to the beginning of each filename.

```bash
# In dir: ~/Projects/Mutant_Barley/longranger_morex_v3/M01-3-3/outs
rename -v "" "M01-3-3_" *

# In dir: ~/Projects/Mutant_Barley/longranger_morex_v3/M20-2-2/outs
rename -v "" "M20-2-2_" *

# In dir: ~/Projects/Mutant_Barley/longranger_morex_v3/M29-2-2/outs
rename -v "" "M29-2-2_" *

# In dir: ~/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/outs
rename -v "" "morex-sample2_" *
```
