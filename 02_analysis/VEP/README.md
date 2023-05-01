# Analyses of 3 mutated barley lines

### VeP

Prepare vcf files if they have not been tabix indexed. VeP only works on bgzipped and tabix indexed VCF files.

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
# Rare vs common SNPs
sbatch Vep-hybrid_rare_snps.sh
sbatch Vep-hybrid_common_snps.sh
# Rare vs common indels
sbatch Vep-hybrid_rare_indels.sh
sbatch Vep-hybrid_common_indels.sh

# Morex-sample2
./Vep_morex-sample2.sh
```

**Note:** When running VeP, don't use `--total_length` flag. Turning on this flag messes up the file format for BAD_Mutations `Vep_to_Subs.py` script.

### Exploration

Compare stop gain SNPs identified in VeP and Annovar.

Pull out stop_gained only form VeP report.

```bash
# Dependencies
module load perl/modules.centos7.5.26.1
module load htslib/1.9
export PATH=$PATH:/panfs/jay/groups/9/morrellp/shared/Software/ensembl-vep-release-108.1

MUT_VEP_REPORT="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/VEP/HC_LC_gff_SNPs_private_all_samples/mut8_and_3mut10xGenomics.SNPs.private.txt"

MUT_PREFIX=$(basename $MUT_VEP_REPORT .txt)
MUT_OUT_DIR=$(dirname $MUT_VEP_REPORT)

filter_vep -i $MUT_VEP_REPORT -o $MUT_OUT_DIR/$MUT_PREFIX.stop_gained.txt -filter "Consequence is stop_gained"

# Get uniq lines
uniq $MUT_OUT_DIR/$MUT_PREFIX.stop_gained.txt > $MUT_OUT_DIR/$MUT_PREFIX.stop_gained.uniq.txt
```
