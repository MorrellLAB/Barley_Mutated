# igv-reports

Check variants using [igv-reports](https://github.com/igvteam/igv-reports).

For small VCF files (less than ~2000 variants), `create_reports` can be run interactively on MSI in under 5-10 min and are documented below. Larger VCFs will be submitted as Slurm job script.

### Dependencies and shared variables

```bash
# Dependencies
module load python3/3.8.3_anaconda2020.07_mamba
conda activate igvreports
```

```bash
# Shared variables
REF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta"

BAM_Morex="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/morex-sample2/outs/morex-sample2_phased_possorted_bam.bam"

BAM1="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M01-3-3/outs/M01-3-3_phased_possorted_bam.bam"
BAM2="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M20-2-2/outs/M20-2-2_phased_possorted_bam.bam"
BAM3="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/M29-2-2/outs/M29-2-2_phased_possorted_bam.bam"

# Other tracks to add
REF_DIFFS_85xONT_INS="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/nanopore_morex/Morex_85x_ont/filtered/parts_pos/morex_85x_ont.noHomRef.geSup5.callable.parts.INS.bed"
```

### Check 10x Genomics large SVs

10x Genomics 3 mutated lines deletions VCF.

```bash
CURR_VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/mut_3_lines_dels.private.callable.noMorexDiffs.vcf.gz"
# Generate report
create_report $CURR_VCF $REF \
    --tracks $CURR_VCF $BAM1 $BAM2 $BAM3 \
    --output ~/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/igv-report_mut_3_lines_dels.private.callable.noMorexDiffs.html
```

10x Genomics 3 mutated lines large SVs VCF inversions.

```bash
CURR_VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/mut_3_lines_large_svs.INV.private.callable.noMorexDiffs.vcf.gz"
# Generate report
create_report $CURR_VCF $REF \
    --tracks $CURR_VCF $BAM3 \
    --output ~/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/igv-report_mut_3_lines_large_svs.INV.private.callable.noMorexDiffs.html
```

10x Genomics 3 mutated lines large SVs VCF deletion.

```bash
CURR_VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/mut_3_lines_large_svs.DEL.private.callable.noMorexDiffs.vcf.gz"
# Generate report
create_report $CURR_VCF $REF \
    --tracks $CURR_VCF $BAM2 \
    --output ~/Projects/Mutant_Barley/longranger_morex_v3/filtered/quality_filtered/igv-report_mut_3_lines_large_svs.DEL.private.callable.noMorexDiffs.html
```
