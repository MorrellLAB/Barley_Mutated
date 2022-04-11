# VCF exploration

Explore VCFs to get a sense of how to filter variants appropriately.

### 10x Genomics, Nanopore, and PacBio data variant visualization

Visualize variants using Samplot.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing
# Prepare BED files for samplot
./prep_SV_regions_Morex_10x_nanopore_pacbio.sh

# Run samplot
sbatch --array=0-1 samplot-Morex_10x_nanopore_pacbio.sh
```

Visualize variants from VCFs.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing
# Reheader and merge VCFs
./reheader_and_merge_vcf-Morex.sh
```
