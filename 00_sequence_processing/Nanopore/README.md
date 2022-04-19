# Nanopore data processing

Alignment and other processing steps for Nanopore data.

---

### Minimap2 and NGMLR to realign complex regions

Originally, we planned to use the Vulcan pipeline (https://gitlab.com/treangenlab/vulcan/-/tree/master) for alignment and realigning complex regions to get better SV calling downstream. But, Vulcan is not very customizable, so we'll follow some of their general steps (Minimap2 and Samtools steps) written as separate bash scripts combined with modifications to their `vulcan` script.

Align with Minimap2 and only keep primary mapping.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/Nanopore
sbatch np_read_mapping-Morex-sample2.job
# Keep primary mapping and add @SQ header lines
sbatch keep_primary_mapping-Morex-sample2.job

# Repeat for 3 mutated barley samples (M01, M20, M29)
sbatch np_read_mapping-M01.job
sbatch np_read_mapping-M20.job
sbatch np_read_mapping-M29.job
```

After modifying the `vulcan` script, transferred the script to MSI. The modification was adding the `--full_sam` and `--full_sam_primary` options to take a custom generated Minimap2 SAM files as input. Modifications are documented in this forked repository under [commit 0dec3094](https://gitlab.com/ChaochihL/vulcan/-/commit/8dce5d4eb75a6044e0fcc00894e22933f56e91c2). The modified `vulcan` script and documentation are available in this forked repository: https://gitlab.com/ChaochihL/vulcan.

```bash
# In dir: ~/GitHub/vulcan
scp vulcan liux1299@mesabi.msi.umn.edu:/panfs/roc/groups/9/morrellp/liux1299/.conda/envs/vulcan_env/bin
```

Run modified vulcan pipeline.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/Nanopore
sbatch vulcan_read_mapping-Morex-sample2_partsRef.job

# Repeat for 3 mutated barley samples (M01, M20, M29)
# except submit as Slurm task array
sbatch --array=0-2 keep_primary_mapping-mut_lines.job
```

The Vulcan pipeline outputs a sorted BAM file.

Add @RG header line to BAM file, many downstream programs require @RG header lines to be present.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/Nanopore
sbatch add_RG_header.sh
```

Run Sniffles2.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/Nanopore
sbatch sniffles-Morex-sample2_partsRef.job
```

Reheader VCF so `SAMPLE` gets replaced with actual sample ID.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/Nanopore
./reheader_vcf-ONT_Morex-sample2.sh
```
