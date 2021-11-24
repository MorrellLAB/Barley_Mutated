# ANNOVAR on mutated barley dataset

Run [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/) to annotate intergenic, exonic, intronic, synonymous, and nonsynonymous SNPs.

**Purpose:** To prepare the nonsynonymous SNPs for predicting deleterious SNPs using BAD_Mutations and other downstream tools.

---

### Files required

- VCF file: `/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated/Filtered/mut_3_lines_filtered_singletons_only.vcf`
- Reference genome file: `/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules_parts.fasta`
- GTF file: 

### File Preparation

Convert GFF to GTF format. Here is a useful comparison of various GFF to GTF conversion tools: https://github.com/NBISweden/GAAS/blob/master/annotation/knowledge/gff_to_gtf.md. We'll go with the recommended AGAT tool, which correctly converts GFF3 to GTF3 format.

First, use Singularity to pull AGAT (available biocontainers are linked on this page https://github.com/NBISweden/AGAT#installation).

```bash
# In dir: ~/Shared/Software
# Get the chosen AGAT container version
singularity pull docker://quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0
# Run the container
# Since we are running this on MSI, we'll need to bind the Shared directory in the container
#   to access the shared directory from within the container
singularity run -B /panfs/roc/groups/9/morrellp/shared:/Shared_container agat_0.8.0--pl5262hdfd78af_0.sif

# Now, we are inside the container and can call on agat the following way
agat_convert_sp_gxf2gxf.pl --help
```

Now, we are inside of the singularity container, let's navigate the the directory containing the GFF files.

```bash
# Morex v2
# Navigate to directory containing GFF files
Singularity> cd /Shared_container/References/Reference_Sequences/Barley/Morex_v2/gene_annotation/

# Convert GFF to GTF
Singularity> agat_convert_sp_gff2gtf.pl --gff Barley_Morex_V2_gene_annotation_PGSB.all.parts.gff3 -o Barley_Morex_V2_gene_annotation_PGSB.all.parts.unsorted.gtf

Singularity> agat_convert_sp_gff2gtf.pl --gff Barley_Morex_V2_gene_annotation_PGSB.HC.parts.gff3 -o Barley_Morex_V2_gene_annotation_PGSB.HC.parts.unsorted.gtf

Singularity> agat_convert_sp_gff2gtf.pl --gff Barley_Morex_V2_gene_annotation_PGSB.LC.gff3 -o Barley_Morex_V2_gene_annotation_PGSB.LC.unsorted.gtf
```

AGAT outputs an unsorted GTF file, but it seems like most downstream tools don't need a sorted GTF. To be safe, we'll sort the GTF file.

```bash
# In dir: /panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/gene_annotation
# Sort GTF files
~/Shared/Software/gff3sort/gff3sort.pl Barley_Morex_V2_gene_annotation_PGSB.all.parts.unsorted.gtf > Barley_Morex_V2_gene_annotation_PGSB.all.parts.sorted.gtf

~/Shared/Software/gff3sort/gff3sort.pl Barley_Morex_V2_gene_annotation_PGSB.HC.parts.unsorted.gtf > Barley_Morex_V2_gene_annotation_PGSB.HC.parts.sorted.gtf

~/Shared/Software/gff3sort/gff3sort.pl Barley_Morex_V2_gene_annotation_PGSB.LC.unsorted.gtf > Barley_Morex_V2_gene_annotation_PGSB.LC.sorted.gtf

# Cleanup unsorted GTF files
rm *.unsorted.gtf
```

*Note:* The reason we are using [gff3sort](https://github.com/billzt/gff3sort) to sort the GTF files instead of UNIX `sort` is because GFF/GTF files have lines with the same chromosomes and start positions. With UNIX `sort`, these lines would be place randomly and sometimes parent features can end up getting place after their children lines. So, we use the [gff3sort](https://github.com/billzt/gff3sort) tool to make sure the GTF is sorted correctly.

### ANNOVAR Steps

Load dependencies for ANNOVAR.

```bash
module load perl/5.26.1
```