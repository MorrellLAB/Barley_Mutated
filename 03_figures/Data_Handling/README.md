# Data handling scripts to parse/prepare results files for plotting

A place to store all data handling/parsing scripts to prepare results files for plotting when they don't have a specific analysis they belong to (e.g., BAD_Mutations results parsing scripts will be stored in `02_analysis/bad_mutations` subdirectory whereas morex-sample2 indel size distribution script will be stored here `03_figures/Data_Handling`).

### morex-sample2 indel sizes

Get size range of indels in morex-sample2.

```bash
module load python3/3.9.3_anaconda2021.11_mamba

# In dir: ~/GitHub/Barley_Mutated/03_figures/Data_Handling
./indel_size_table.py \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/diffs_from_ref/morex-sample2_phased_variants-indels.callable.biallelic.vcf.gz \
    /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/diffs_from_ref \
    "morex-sample2"
```

### Mutated samples indels

Select indels 1-3 bp in size.

```bash
module load python3/3.8.3_anaconda2020.07_mamba

# In dir: ~/GitHub/Barley_Mutated/03_figures/Data_Handling
./select_indels_by_size.py /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.INDELs.private.vcf.gz 1 3 > ~/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.INDELs.private.1-3bp.vcf

# Per sample counts of 1-3 bp indels
./count_indels_per_individual.py /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.INDELs.private.1-3bp.vcf ~/Projects/Mutant_Barley/tables mut8_and_3mut10xGenomics.INDELs.private.1-3bp
```
