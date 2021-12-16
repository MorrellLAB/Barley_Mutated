# ANNOVAR on mutated barley dataset

Run [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/) to annotate intergenic, exonic, intronic, synonymous, and nonsynonymous SNPs.

**Purpose:** To prepare the nonsynonymous SNPs for predicting deleterious SNPs using BAD_Mutations and other downstream tools.

---

### Files required

- VCF file: `/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated/Filtered/mut_3_lines_filtered_singletons_only.vcf`
- Reference genome file: `/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules_parts.fasta`
- GFF3 file(s): `/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/gene_annotation/Barley_Morex_V2_gene_annotation_PGSB.HC.parts.gff3`

### ANNOVAR Steps

We will mostly follow the guide here: https://annovar.openbioinformatics.org/en/latest/user-guide/gene/.

Load dependencies for ANNOVAR.

```bash
module load perl/5.26.1

# Export paths to directories containing annovar scripts so scripts are calleble from anywhere without specifying the path
export PATH=${PATH}:/panfs/roc/groups/9/morrellp/shared/Software/annovar
export PATH=${PATH}:/panfs/roc/groups/9/morrellp/shared/Software/annovar_conversion_tools
```

Generate the genePred file from the GFF3 file.

```bash
# In dir: ~/Projects/Mutant_Barley/results/Annovar/HC
gff3ToGenePred ~/Shared/References/Reference_Sequences/Barley/Morex_v2/gene_annotation/Barley_Morex_V2_gene_annotation_PGSB.HC.parts.gff3 HV_Morex_v2_HC_refGene0.txt
# Sort genePred file
sort -V -k1,1 HV_Morex_v2_HC_refGene0.txt > HV_Morex_v2_HC_refGene0_sorted.txt
```

Add a random first field to the genePred file using nl (this add the line number to each line).

```bash
nl HV_Morex_v2_HC_refGene0_sorted.txt > HV_Morex_v2_HC_refGene.txt

# Cleanup intermediate files
rm HV_Morex_v2_HC_refGene0.txt
rm HV_Morex_v2_HC_refGene0_sorted.txt
```

Build the database: generate a transcript FASTA file.

```bash
# Define some variables for filepaths
ref_fasta="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules_parts.fasta"
genePred_file="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/Annovar/HC/HV_Morex_v2_HC_refGene.txt"
out_prefix="HV_Morex_v2_HC_refGene"

# Run annovar script
retrieve_seq_from_fasta.pl --format refGene --seqfile ${ref_fasta} ${genePred_file} --outfile ${out_prefix}Mrna.fa
```

Prepare VCF file for Annovar. Convert VCF to Annovar's input format using Annovar's `convert2annovar.pl` script. We'll need to submit this as a job to MSI.

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/Annovar
sbatch vcf_to_annovar_input.job
```

Annotation with Annovar using the script `annotate_variation.pl`. This will also need to be submitted as a job to MSI.

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/Annovar
sbatch annotate_with_annovar.job
```

### Steps to prepare for BAD_Mutations run

Separate insertions and deletions into a separate file from SNPs, SNPs will be run separately through BAD_Mutations.

```bash
# In dir: ~/Projects/Mutant_Barley/results/Annovar/HC
# Get a sense of the unqiue types of annotations in the file
cut -f 2 mut_3_lines_filtered_singletons_only_annovar_input.txt.exonic_variant_function | sort -uV
# Output
frameshift deletion
frameshift insertion
nonframeshift deletion
nonsynonymous SNV
stopgain
stoploss
synonymous SNV

cut -f 1 mut_3_lines_filtered_singletons_only_annovar_input.txt.variant_function | sort -uV
# Output
downstream
exonic
intergenic
intronic
splicing
upstream
upstream;downstream

# Check that we're only pulling insertions and deletions
grep -w 'deletion\|insertion' mut_3_lines_filtered_singletons_only_annovar_input.txt.exonic_variant_function | cut -f 2 | sort -uV
# Output
frameshift deletion
frameshift insertion
nonframeshift deletion

# Check that we're only pulling everything that's not an indel
grep -vw 'deletion\|insertion' mut_3_lines_filtered_singletons_only_annovar_input.txt.exonic_variant_function | cut -f 2 | sort -uV
# Output
nonsynonymous SNV
stopgain
stoploss
synonymous SNV

# Save to file
# Indels
grep -w 'deletion\|insertion' mut_3_lines_filtered_singletons_only_annovar_input.txt.exonic_variant_function > mut_3_lines_filtered_singletons_only_annovar_input.txt.indels.exonic_variant_function
# SNVs
grep -vw 'deletion\|insertion' mut_3_lines_filtered_singletons_only_annovar_input.txt.exonic_variant_function > mut_3_lines_filtered_singletons_only_annovar_input.txt.SNVs.exonic_variant_function
```

We'll also generate some toy files of the ANNOVAR output and put it in a subdirectory called `toy_datasets` so it's easier to track file format changes in the future and modify the `ANNOVAR_To_Effects.py` script accordingly.

```bash
# In dir: ~/Projects/Mutant_Barley/results/Annovar/HC
# Generate toy files
head -n 20 mut_3_lines_filtered_singletons_only_annovar_input.txt.SNVs.exonic_variant_function > toy.SNVs.exonic_variant_function

head -n 20 mut_3_lines_filtered_singletons_only_annovar_input.txt.indels.exonic_variant_function > toy.indels.exonic_variant_function

# From the variant_function file, we'll want to make sure we include all variants present in the indels and SNVs file above
head -n 60 mut_3_lines_filtered_singletons_only_annovar_input.txt.variant_function > toy.variant_function
```

Determine whether and where pre-mature stop codons could be located. See documentation for details: https://annovar.openbioinformatics.org/en/latest/misc/accessory/#coding_change-infer-mutated-protein-sequence

```bash
# Define variables
exon_var_fn_file="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/Annovar/HC/mut_3_lines_filtered_singletons_only_annovar_input.txt.exonic_variant_function"
ref_gene_file="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/Annovar/HC/HV_Morex_v2_HC_refGene.txt"
ref_gene_mrna_fa_file="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/Annovar/HC/HV_Morex_v2_HC_refGeneMrna.fa"
# Run script
coding_change.pl "${exon_var_fn_file}" "${ref_gene_file}" "${ref_gene_mrna_fa_file}"
```

Make a unified table using the script: `ANNOVAR_To_Effects.py` (Note: this script needs to be modified from Tom and Li's previous versions to handle insertions and deletions correctly). This unified table will be used in BAD_Mutations.

Make a unified table of just the SNVs first.

```bash
# In dir: ~/Projects/Mutant_Barley/results/Annovar/HC
module load python3/3.8.3_anaconda2020.07_mamba

# Define some variables
var_fun="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/Annovar/HC/mut_3_lines_filtered_singletons_only_annovar_input.txt.variant_function"
exon_var_fun="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/Annovar/HC/mut_3_lines_filtered_singletons_only_annovar_input.txt.SNVs.exonic_variant_function"
out_prefix="mut_3_lines_filtered_singletons_only_annovar_SNVs"

~/GitHub/Barley_Mutated/02_analysis/Annovar/ANNOVAR_To_Effects.py ${var_fun} ${exon_var_fun} > ${out_prefix}_unified.table
```
