# smoove

Smoove is intented to assist with small/mid-sized SV calling using short read data. We'll run [smoove](https://github.com/brentp/smoove) on the 10x Genomics three mutated lines and morex-sample2, 8 WGS mutated lines, and hybrid barley landraces to detect small and mid-sized SVs. For the 10x Genomics dataset, we can intersect the SVs called to see which SVs are called by both pipelines. We'll still rely on GATK for SNP calls from the short read data.

### Methods

Run step 1 of population-level calling (large cohorts and parallelize across samples utilizing Slurm job arrays.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/smoove
sbatch --array=0-3 run_smoove-step1.sh
# Some samples (array indices) may timeout or fail, figure out which ones
#   to re-run using the get_re-run_array_indices.sh script available here:
#   https://github.com/MorrellLAB/Barley_Mutated/blob/master/02_analysis/bad_mutations/get_re-run_array_indices.sh
get_re-run_array_indices.sh 146875535
0-1
# Increase walltime and re-run timeout indices
sbatch --array=0-1 run_smoove-step1.sh
```
