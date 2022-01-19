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

After modifying the scripts, use `sbatch` to submit each script.
