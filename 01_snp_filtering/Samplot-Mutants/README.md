# VCF exploration

Final check of filtered 10x Genomics mutated lines.

### Score post-filtering SVs

Generate Samplot images.

```bash
# In dir: ~/GitHub/Barley_Mutated/01_snp_filtering/Samplot-Mutants
sbatch samplot_vcf-mut_10x.sh
```

After filtering (to reduce total number of SVs), generate Samplot images for SV-plaudit. We'll score the 155 DELs from the file `/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/de_novo_larger_svs/mut_3_lines_dels_merged.callable.noRefDiffs.private.vcf.gz` to finalize filtering this file.

Run SV-Plaudit following PlotCritic setup instructions in Github repo: https://github.com/jbelyeu/SV-plaudit.

```bash
# Load dependencies
module load python3/3.8.3_anaconda2020.07_mamba

# Create a PlotCritic website
# Ran the following substituting our own fields
python /panfs/jay/groups/9/morrellp/liux1299/Software/SV-plaudit/PlotCritic/project_setup.py \
    -p mut_10x_Genomics \
    -e example_user@umn.edu \
    -a [ACCESS_KEY] -s [SECRET_ACCESS_ID] \
    -q "Does evidence in the sample support the variant called?" \
    -A "s":"Supports" "n":"Does not support" "d":"De novo" -r \
    -R "chrom" "start" "end" "sv_type" "titles" "bams" \
    -S "chrom" "start" "end" "sv_type"

# Upload images to PlotCritic website
python /panfs/jay/groups/9/morrellp/liux1299/Software/SV-plaudit/PlotCritic/upload.py \
    -d /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/samplot-mut_10x/DEL \
    -c /panfs/jay/groups/9/morrellp/liux1299/Software/SV-plaudit/PlotCritic/config.json
```

Score images, then retrieve scores and pull out supported SVs. Download both `summary_report.tsv` and `raw_report.tsv`, then rename files to include sample name. *Note:* Tried the command line retrieval but that resulted in an empty file, so we just downloaded it from the PlotCritic site where we scored the variants.

Couldn't get SV-plaudit's `annotate.py` script to work (I think there are multiple one-off errors with the field that parts of the script is pulling from). So, we'll use a workaround instead.

```bash
# In dir: ~/Projects/Mutant_Barley/samplot-mut_10x
# Create file containing only supports positions
awk '$5 == "100.0" { print $0 }' summary_report-mut_10x_dels.tsv | cut -f 1-3 | sort -k1,1 -k2,2n > mut_10x_dels.supports.bed

# Use bedtools
module load bedtools/2.29.2
# Prepare header
zgrep "#" /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/de_novo_larger_svs/mut_3_lines_dels_merged.callable.noRefDiffs.private.vcf.gz > mut_10x_dels.callable.noRefDiffs.private.supports.vcf
# Intersect VCF with supports bed file
bedtools intersect -wa -a /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/de_novo_larger_svs/mut_3_lines_dels_merged.callable.noRefDiffs.private.vcf.gz -b mut_10x_dels.supports.bed | sort -uV -k1,1 -k2,2n >> mut_10x_dels.callable.noRefDiffs.private.supports.vcf
```

VCF file including only variants that were scored as "supports":

```bash
/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/samplot-mut_10x/mut_10x_dels.callable.noRefDiffs.private.supports.vcf
```
