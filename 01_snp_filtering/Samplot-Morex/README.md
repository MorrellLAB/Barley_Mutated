# VCF exploration

Explore VCFs to get a sense of how to filter variants appropriately.

### 10x Genomics, Nanopore, and PacBio data variant visualization

Visualize variants using Samplot.

We'll generate images from the VCF, this requires all samples be in a single VCF file.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/Samplot-Morex
merge_vcfs-morex_10x_ont_pacbio.sh
```

Samplot images of merged VCF 10x ONT PacBio.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/Samplot-Morex
sbatch samplot_vcf-Morex_10x_ont_pacbio.sh
```

Samplot images of smoove output VCF that includes Morex 10x, 3 mutated lines 10x, 8 WGS mutated lines, and WGS hybrid barley lines.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/Samplot-Morex
sbatch --array=0-4 samplot_vcf-smoove_Morex10x-wgs_mut-hybrid.sh
```

Next, generate important regions from the VCF containing ONT and PacBio data relative to Morex v3.

```bash
# In dir: ~/Projects/Mutant_Barley/exploration_morex_v3
module load python3/3.8.3_anaconda2020.07_mamba
~/GitHub/Barley_Mutated/01_snp_filtering/vcf_sniffles_long_read_to_bed.py morex_ont_pacbio_noHomRef.vcf > morex_ont_pacbio_noHomRef_noBND.bed
```

Get BED files for each variant type individually, this is so windows are only for the variant type of interest and not any other variant type.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/Samplot-Morex
./prep_SV_regions_Morex_ont_pacbio.sh
./prep_SV_regions_filt_Morex_10x_ont_pacbio.sh
```

Convert TE annotations to BED format.

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v3/repeat_annotation
module load bedops_ML/2.4.38
gff2bed < TEanno-v1.0__200416_MorexV3_pseudomolecules.sorted.gff > TEanno-v1.0__200416_MorexV3_pseudomolecules.sorted.bed
cut -f 1-3 TEanno-v1.0__200416_MorexV3_pseudomolecules.sorted.bed > TEanno-v1.0__200416_MorexV3_pseudomolecules.sorted.minimal.bed
# Convert to parts positions
~/GitHub/File_Conversions/Barley_Pseudomolecules_to_Parts.py --gff TEanno-v1.0__200416_MorexV3_pseudomolecules.sorted.gff morex_v3 > TEanno-v1.0__200416_MorexV3_pseudomolecules.sorted.parts.gff
gff2bed < TEanno-v1.0__200416_MorexV3_pseudomolecules.sorted.parts.gff > TEanno-v1.0__200416_MorexV3_pseudomolecules.sorted.parts.bed
cut -f 1-3 TEanno-v1.0__200416_MorexV3_pseudomolecules.sorted.parts.bed > TEanno-v1.0__200416_MorexV3_pseudomolecules.sorted.parts.minimal.bed
bgzip -c TEanno-v1.0__200416_MorexV3_pseudomolecules.sorted.parts.minimal.bed > TEanno-v1.0__200416_MorexV3_pseudomolecules.sorted.parts.minimal.bed.gz
tabix -p bed TEanno-v1.0__200416_MorexV3_pseudomolecules.sorted.parts.minimal.bed.gz
```

Run samplot.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/Samplot-Morex
# Run samplot
sbatch --array=0-2 samplot_plot-Morex_ont_pacbio.sh
sbatch --array=0-2 samplot_plot-Morex_10x_ont_pacbio.sh
```

Visualize variants from VCFs.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing
# Reheader and merge VCFs
./reheader_and_merge_vcf-Morex.sh
```
