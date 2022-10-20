# smoove

Smoove is intented to assist with small/mid-sized SV calling using short read data. We'll run [smoove](https://github.com/brentp/smoove) on the 10x Genomics three mutated lines and morex-sample2, 8 WGS mutated lines, and hybrid barley landraces to detect small and mid-sized SVs. For the 10x Genomics dataset, we can intersect the SVs called to see which SVs are called by both pipelines. We'll still rely on GATK for SNP calls from the short read data.

### Methods

**Prepare exclusion regions:** Mask problematic regions using sample-specific genome exclude files. Regions with depth much larger than the norm are problematic as well as gaps in the reference genome, these can produce low confidence or misleading SV calls.

Calculate the depth of each BAM file with mosdepth.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/smoove
sbatch --array=0-38 run_mosdepth.sh
```

Extract high depth regions.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/smoove
sbatch run_extract_high_depth_regions.sh
```

For each sample specific high depth region, concat and merge with gap regions.

```bash
sbatch --array=0-38 create_exclude_files.sh
```

Run step 1 of population-level calling (large cohorts and parallelize across samples utilizing Slurm job arrays.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/smoove
sbatch --array=0-38 run_smoove-step1.sh
# Some samples (array indices) may timeout or fail, figure out which ones
#   to re-run using the get_re-run_array_indices.sh script available here:
#   https://github.com/MorrellLAB/Barley_Mutated/blob/master/02_analysis/bad_mutations/get_re-run_array_indices.sh
# Except the ones below, all other samples completed within 36 hours
get_re-run_array_indices.sh 147219025
    3
# Increase walltime and re-run timeout indices
sbatch --array=3 run_smoove-step1.sh
```

Run steps 2-5 of population-level calling together utilizing GNU parallel where appropriate.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/smoove
sbatch run_smoove-step2-5.sh
```

Separate smoove called samples into 10x Genomics mutated lines, 8 mutated wgs lines, and the 13 hybrid barley lines relevant to this study.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/smoove
sbatch select_samples_smoove.sh
```
