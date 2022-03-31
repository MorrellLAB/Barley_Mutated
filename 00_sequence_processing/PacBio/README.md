# PacBio data processing

We used PacBio data of Morex published in Mascher et al. 2021. Here is a summary of the SRA numbers for the data we downloaded.

| Description | BioProject | ERS numbers | ERX numbers | ERR numbers |
|-------------|------------|-------------|-------------|-------------|
| WGS Morex V3 CCS raw reads | [PRJEB40587](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB40587/) | ERS5135541-ERS5135545 | ERX4582061-ERX4582065 | ERR4659245-ERR4659249 |

### Download data

Download SRA files.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/PacBio
sbatch sra_download.sh
```

Convert SRA to FASTQ format.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/PacBio
sbatch sra_to_fastq.sh
```

### Sample QC

Run LongQC to do quality assessment of FASTQ files.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/PacBio
sbatch quality_assessment_raw_fastq.sh 
```

### Sample Processing

Adapter filtering with [HiFiAdapterFilt](https://github.com/sheinasim/HiFiAdapterFilt) tool.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/PacBio
sbatch adapter_filtering_pacbio.sh
```

Concatenate FASTQ files into single file before aligning.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/PacBio
sbatch concat_fastq_Morex_pacbio.job
```

Align using Minimap2.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/PacBio
sbatch minimap2_read_mapping-Morex_pacbio.job

# Keep primary mapping and add @SQ header lines
sbatch keep_primary_mapping-Morex_pacbio.job
```

Run modified vulcan pipeline. The modification was adding the `--ful_sam` and `--full_sam_primary` options to take a custom generated Minimap2 SAM files as input. Modifications are documented in this forked repository under [commit 0dec3094](https://gitlab.com/ChaochihL/vulcan/-/commit/8dce5d4eb75a6044e0fcc00894e22933f56e91c2). The modified `vulcan` script and documentation are available in this forked repository: https://gitlab.com/ChaochihL/vulcan.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/PacBio
sbatch vulcan_read_mapping-Morex_pacbio.job
```

The Vulcan pipeline outputs a sorted BAM file.

Run Sniffles2.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/PacBio
sbatch sniffles-Morex_pacbio.sh
```

Add @RG header line to BAM file, many downstream programs require @RG header lines to be present.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/PacBio
sbatch add_RG_header.sh
```
