# ANNOVAR on mutated barley dataset

Run [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/) to annotate intergenic, exonic, intronic, synonymous, and nonsynonymous SNPs.

**Purpose:** To prepare the nonsynonymous SNPs for predicting deleterious SNPs using BAD_Mutations and other downstream tools.

---

### Files required

- VCF file: `/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.SNPs.private.vcf.gz`
- Reference genome file: `/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta`
- GFF3 file(s):
    - HC: `/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/gene_annotation/Hv_Morex.pgsb.Jul2020.HC.sorted.parts.gff3`
    - All: `/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/gene_annotation/Hv_Morex.pgsb.Jul2020.sorted.parts.gff3`

### ANNOVAR Steps

We will mostly follow the guide here: https://annovar.openbioinformatics.org/en/latest/user-guide/gene/. We'll run it for both high confidence gene models and all gene models.

Load dependencies for ANNOVAR.

```bash
module load perl/5.26.1
# Export paths to directories containing annovar scripts so scripts are calleble from anywhere without specifying the path
export PATH=${PATH}:/panfs/jay/groups/9/morrellp/shared/Software/annovar
export PATH=${PATH}:/panfs/jay/groups/9/morrellp/shared/Software/annovar_conversion_tools
```

---

#### Run Annovar on HC only and all gene models

Filepaths

```bash
GFF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/gene_annotation/Hv_Morex.pgsb.Jul2020.sorted.parts.gff3"
GFF_HC="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/gene_annotation/Hv_Morex.pgsb.Jul2020.HC.sorted.parts.gff3"

OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/Annovar/all"
OUT_DIR_HC="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/Annovar/HC"

REF_FASTA="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta"

MUT_VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.SNPs.private.vcf.gz"
MUT_VCF_INDELs="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.INDELs.private.vcf.gz"
```

Generate the genePred file from the GFF3 file.

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/Annovar
# All (HC and LC) genes
./gff3_to_gene_pred.sh ${GFF} ${OUT_DIR} hv_morex_v3_all
# HC only
./gff3_to_gene_pred.sh ${GFF_HC} ${OUT_DIR_HC} hv_morex_v3_HC
```

Build the database: generate a transcript FASTA file.

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/Annovar
# Run annovar script
# All (HC and LC) genes
./build_db.sh ${OUT_DIR}/hv_morex_v3_all_refGene.txt ${REF_FASTA} ${OUT_DIR} hv_morex_v3_all
# HC only
./build_db.sh ${OUT_DIR_HC}/hv_morex_v3_HC_refGene.txt ${REF_FASTA} ${OUT_DIR_HC} hv_morex_v3_HC
```

Prepare VCF file for Annovar. Convert VCF to Annovar's input format using Annovar's `convert2annovar.pl` script.

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/Annovar
# SNPs
./vcf_to_annovar_input.sh ${MUT_VCF} ${OUT_DIR}
./vcf_to_annovar_input.sh ${MUT_VCF} ${OUT_DIR_HC}

# indels
./vcf_to_annovar_input.sh ${MUT_VCF_INDELs} ${OUT_DIR}
./vcf_to_annovar_input.sh ${MUT_VCF_INDELs} ${OUT_DIR_HC}
```

Annotation with Annovar using the script `annotate_variation.pl`.

```bash
# In dir: ~/GitHub/Barley_Mutated/02_analysis/Annovar
# SNPs
./annotate_with_annovar.sh "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/Annovar/all/mut8_and_3mut10xGenomics.SNPs.private_annovar_input.txt" ${OUT_DIR} "hv_morex_v3_all"

./annotate_with_annovar.sh "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/Annovar/HC/mut8_and_3mut10xGenomics.SNPs.private_annovar_input.txt" ${OUT_DIR_HC} "hv_morex_v3_HC"

# indels
./annotate_with_annovar.sh "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/Annovar/all/mut8_and_3mut10xGenomics.INDELs.private_annovar_input.txt" ${OUT_DIR} "hv_morex_v3_all"

./annotate_with_annovar.sh "/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/Annovar/HC/mut8_and_3mut10xGenomics.INDELs.private_annovar_input.txt" ${OUT_DIR_HC} "hv_morex_v3_HC"
```

---

### Steps to prepare for BAD_Mutations run

**No longer used, ended up using VeP output instead but will leave notes here for future reference.**

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

# Check that we're only pulling everything that's not an indel or a stopgain
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
# stopgain only
grep -w 'stopgain' mut_3_lines_filtered_singletons_only_annovar_input.txt.exonic_variant_function > mut_3_lines_filtered_singletons_only_annovar_input.txt.stopgain.exonic_variant_function
# Create a list of gene names so we can crosscheck in case we run into errors during the BAD_Mutations align step
cut -f 3 mut_3_lines_filtered_singletons_only_annovar_input.txt.stopgain.exonic_variant_function | tr ':' '\t' | cut -f 2 | sort -uV > stopgain_gene_list.txt
```

We'll also generate some toy files of the ANNOVAR output and put it in a subdirectory called `toy_datasets` so it's easier to track file format changes in the future and modify the `ANNOVAR_To_Effects.py` script accordingly.

```bash
# In dir: ~/Projects/Mutant_Barley/results/Annovar/HC
# Generate toy files
head -n 20 mut_3_lines_filtered_singletons_only_annovar_input.txt.SNVs.exonic_variant_function > toy_HC.SNVs.exonic_variant_function

head -n 20 mut_3_lines_filtered_singletons_only_annovar_input.txt.indels.exonic_variant_function > toy_HC.indels.exonic_variant_function

# From the variant_function file, we'll want to make sure we include all variants present in the indels and SNVs file above
head -n 60 mut_3_lines_filtered_singletons_only_annovar_input.txt.variant_function > toy_HC.variant_function
```

Make a unified table using the script: `ANNOVAR_To_Effects.py` (Note: this script needs to be modified from Tom and Li's previous versions to handle insertions and deletions correctly). This unified table will be used in BAD_Mutations.

Make a unified table of just the SNVs first.

```bash
# In dir: ~/Projects/Mutant_Barley/results/Annovar/HC
module load python3/3.8.3_anaconda2020.07_mamba

# Define some variables
var_fun="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/Annovar/HC/mut_3_lines_filtered_singletons_only_annovar_input.txt.variant_function"
exon_var_fun="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/Annovar/HC/mut_3_lines_filtered_singletons_only_annovar_input.txt.SNVs.exonic_variant_function"
out_prefix="mut_3_lines_filtered_singletons_only_annovar_SNVs"

~/GitHub/Barley_Mutated/02_analysis/Annovar/ANNOVAR_To_Effects.py ${var_fun} ${exon_var_fun} > ${out_prefix}_unified.table
```
