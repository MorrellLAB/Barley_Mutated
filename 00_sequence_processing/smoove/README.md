# smoove

Smoove is intented to assist with small/mid-sized SV calling using short read data. We'll run [smoove](https://github.com/brentp/smoove) on the 10x Genomics three mutated lines and morex-sample2, 8 WGS mutated lines, and hybrid barley landraces to detect small and mid-sized SVs. For the 10x Genomics dataset, we can intersect the SVs called to see which SVs are called by both pipelines. We'll still rely on GATK for SNP calls from the short read data.

### Methods

Run step 1 of population-level calling (large cohorts and parallelize across samples utilizing Slurm job arrays.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/smoove
sbatch --array=0-3 run_smoove-step1.sh
```
