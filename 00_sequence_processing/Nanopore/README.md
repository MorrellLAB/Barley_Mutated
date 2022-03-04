# Nanopore data processing

Alignment and other processing steps for Nanopore data.

---

### Minimap2 and NGMLR to realign complex regions

Originally, we planned to use the Vulcan pipeline (https://gitlab.com/treangenlab/vulcan/-/tree/master) for alignment and realigning complex regions to get better SV calling downstream. But, Vulcan is not very customizable, so we'll follow their general steps but written as separate bash scripts with custom modifications.

Align with Minimap2 and only keep primary mapping.

```bash
# In dir: ~/GitHub/Barley_Mutated/00_sequence_processing/Nanopore
sbatch np_read_mapping-Morex-sample2.job

# Keep primary mapping and add @SQ header lines
sbatch keep_primary_mapping-Morex-sample2.job
```
