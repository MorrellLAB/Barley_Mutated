# Analyses of 3 mutated barley lines

### VeP

Prepare vcf files. VeP only works on bgzipped and tabix indexed VCF files.

```bash
# In dir: ~/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated/Filtered
module load htslib/1.9

# bgzip and index
bgzip -c mut_3_lines_filtered_singletons_only.vcf > mut_3_lines_filtered_singletons_only.vcf.gz
tabix -p vcf mut_3_lines_filtered_singletons_only.vcf.gz

bgzip -c mut_3_lines_filtered_hom_singletons_only.vcf > mut_3_lines_filtered_hom_singletons_only.vcf.gz
tabix -p vcf mut_3_lines_filtered_hom_singletons_only.vcf.gz
```

Run VeP with scripts:

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/VEP
# All samples combined
./Vep_singletons-mut_snps.sh
./Vep_singletons-mut_indels.sh
# Vep for each sample individually
./Vep_singletons_by_sample-mut_snps.sh
./Vep_singletons_by_sample-mut_indels.sh

# Hybrid VCF is much larger, will need to submit as Slurm job
sbatch Vep-hybrid_snps.sh
sbatch Vep-hybrid_indels.sh

./Vep_morex-sample2.sh
```

**Note:** When running VeP, don't use `--total_length` flag. Turning on this flag messes up the file format for BAD_Mutations `Vep_to_Subs.py` script.

#### Preparing VeP files for bad_mutations analysis

Pull out missense variants and stop lost variants for BAD_Mutations analysis. BAD_Mutations likely can only annotate missense variants, but might take stop lost variants too, so we will try.

```bash
# In dir: ~/Projects/Mutant_Barley/results/VEP/HC_parts_gff
mkdir for_bad_mut

vep_report=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/HC_parts_gff/mut_3_lines_filtered_singletons_only.txt

prefix=$(basename ${vep_report} .txt)
grep "#" ${vep_report} > for_bad_mut/${prefix}_missense_and_stoplost.txt
# Store each type of consequence in sepearate temp files
grep "missense" ${vep_report} >> tmp_${prefix}_missense.txt
grep "stop_lost" ${vep_report} >> tmp_${prefix}_stoplost.txt
# Combine temp files
cat tmp_${prefix}_missense.txt tmp_${prefix}_stoplost.txt | sort -Vu -k2,2 >> for_bad_mut/${prefix}_missense_and_stoplost.txt
# Vep_to_Subs.py only works with gzipped files
gzip for_bad_mut/${prefix}_missense_and_stoplost.txt
# Cleanup
rm tmp_${prefix}_missense.txt
rm tmp_${prefix}_stoplost.txt
```

**Note:** M01-3-3 doesn't have missense variants or stop lost variants.
