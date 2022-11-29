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

---

### Score post-filtering SVs

#### 10x Genomics morex-sample2

After filtering (to reduce total number of SVs), generate Samplot images for SV-plaudit. We'll score the 208 DELs from the file `/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/morex-sample2_dels.10xCustomFilt.noBND.noRepeatOverlap.noRefNs.vcf.gz` to finalize filtering this file.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/Samplot-Morex
sbatch samplot_vcf-Morex_10x.sh
```

Run SV-Plaudit following PlotCritic setup instructions in Github repo: https://github.com/jbelyeu/SV-plaudit.

```bash
# Load dependencies
module load python3/3.8.3_anaconda2020.07_mamba

# Create a PlotCritic website
# Ran the following substituting our own fields
python /panfs/jay/groups/9/morrellp/liux1299/Software/SV-plaudit/PlotCritic/project_setup.py \
    -p morex_10x_Genomics \
    -e example_user@umn.edu \
    -a [ACCESS_KEY] -s [SECRET_ACCESS_ID] \
    -q "Does evidence in the sample support the variant called?" \
    -A "s":"Supports" "n":"Does not support" "d":"De novo" -r \
    -R "chrom" "start" "end" "sv_type" "titles" "bams" \
    -S "chrom" "start" "end" "sv_type"

# Upload images to PlotCritic website
python /panfs/jay/groups/9/morrellp/liux1299/Software/SV-plaudit/PlotCritic/upload.py \
    -d /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/samplot-morex_10x/DEL \
    -c /panfs/jay/groups/9/morrellp/liux1299/Software/SV-plaudit/PlotCritic/config.json
```

Score images, then retrieve scores and pull out supported SVs. Download both `summary_report.tsv` and `raw_report.tsv`, then rename files to include sample name. *Note:* Tried the command line retrieval but that resulted in an empty file, so we just downloaded it from the PlotCritic site where we scored the variants.

Couldn't get SV-plaudit's `annotate.py` script to work (I think there are multiple one-off errors with the field that parts of the script is pulling from). So, we'll use a workaround instead.

```bash
# In dir: ~/Projects/Mutant_Barley/samplot-morex_10x
# Create file containing only supports positions
awk '$5 == "100.0" { print $0 }' morex_sample2_DEL_summary_report.tsv | cut -f 1-3 | sort -k1,1 -k2,2n > morex_sample2_DEL.supports.bed

# Use bedtools
module load bedtools/2.29.2
# Prepare header
zgrep "#" /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/morex-sample2_dels.10xCustomFilt.noBND.noRepeatOverlap.noRefNs.vcf.gz > ~/Projects/Mutant_Barley/samplot-morex_10x/morex-sample2_dels.10xCustomFilt.noBND.noRepeatOverlap.noRefNs.supports.vcf
# Intersect VCF with supports bed file
bedtools intersect -wa -a /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/morex-sample2_dels.10xCustomFilt.noBND.noRepeatOverlap.noRefNs.vcf.gz -b morex_sample2_DEL.supports.bed | sort -uV -k1,1 -k2,2n >> ~/Projects/Mutant_Barley/samplot-morex_10x/morex-sample2_dels.10xCustomFilt.noBND.noRepeatOverlap.noRefNs.supports.vcf
```

VCF file including only variants that were scored as "supports":

```bash
/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/samplot-morex_10x/morex-sample2_dels.10xCustomFilt.noBND.noRepeatOverlap.noRefNs.supports.vcf
```

#### morex-sample2 ONT

Repeat similar steps above for morex-sample2 ONT data after filtering. After filtering (to reduce total number of SVs), generate Samplot images for SV-plaudit. Browse through images quickly and decide if scoring in SV-plaudit is necessary.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/Samplot-Morex
sbatch samplot_plot-Morex_ont.sh
```

SV-plaudit scoring wasn't necessary here.

#### morex 85x ONT

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/Samplot-Morex
sbatch samplot_plot-Morex_85x_ont.sh
```

SV-plaudit scoring wasn't necessary here.

#### morex PacBio

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/Samplot-Morex
sbatch samplot_plot-Morex_pacbio.sh
```

SV-plaudit scoring wasn't necessary here.

---

## Final checks

Did a final check of 5-10 variants in IGV through MSI's Open On Demand after filtering.
