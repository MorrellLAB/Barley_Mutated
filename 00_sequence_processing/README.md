# Sequence Processing

Mutated barley were aligned to barley parts reference with 10x Genomics software, `longranger`, by John Garbe.
Parts reference filepath: `/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex/barley_RefSeq_v1.0/barley_pseudomolecules_parts.fa`

Summarized coverage using `Coverage_Mapping` from the [`sequence_handling` pipeline](https://github.com/MorrellLAB/sequence_handling).

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
