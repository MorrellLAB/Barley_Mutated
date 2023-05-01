# SnpEff

Run SnpEff as a check/comparison to premature stops identifed by VeP and Annovar.

Run SnpEff on mutated private SNPs.

Convert parts to pseudo positions.

```bash
# ~/Projects/Mutant_Barley/de_novo_vcfs
~/GitHub/File_Conversions/Barley_Parts_to_Pseudomolecules.py --vcf mut8_and_3mut10xGenomics.SNPs.private.vcf.gz morex_v3 > mut8_and_3mut10xGenomics.SNPs.private.pseudo_pos.vcf
```

Use Morex v3 Hordeum_vulgare database from SnpEff.

```bash
# In dir: ~/Software/snpEff
# Dependencies
module load java/openjdk-17.0.2
SNPEFF_JAR="/panfs/jay/groups/9/morrellp/liux1299/Software/snpEff/snpEff.jar"

VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.SNPs.private.pseudo_pos.vcf"
PREFIX="mut8_and_3mut10xGenomics.SNPs.private.pseudo_pos"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/snpeff"

# Annotate VCF
java -jar ${SNPEFF_JAR} -v -stats ${OUT_DIR}/${PREFIX}.html Hordeum_vulgare ${VCF} > ${OUT_DIR}/${PREFIX}.ann.vcf
```





GFF3 to GTF conversion: https://agat.readthedocs.io/en/latest/gff_to_gtf.html

```bash
module load singularity/current
# get the chosen AGAT container version
singularity pull docker://quay.io/biocontainers/agat:1.0.0--pl5321hdfd78af_0
# run the container
singularity run agat_1.0.0--pl5321hdfd78af_0.sif

# Once container is loaded, can't access shared space filepath
cp /panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/gene_annotation/Hv_Morex.pgsb.Jul2020.sorted.parts.gff3 ~

# Convert GFF3 to GTF based on Ensembl's format specification
# https://useast.ensembl.org/info/website/upload/gff.html
agat_convert_sp_gff2gtf.pl --gtf_version 2.2 --gff /home/morrellp/liux1299/Hv_Morex.pgsb.Jul2020.sorted.parts.gff3 --gtf ~/Software/snpEff/data/Hordeum_vulgare_v3_parts/genes.gtf
```

Copy fasta to snpEff data directory.

```bash
cp /panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/gene_annotation/Hv_Morex.pgsb.Jul2020.cds.fa ~/Software/snpEff/data/Hordeum_vulgare_v3_parts/sequences.fa
```

```bash
# Dependencies
module load java/openjdk-17.0.2
SNPEFF_JAR="/panfs/jay/groups/9/morrellp/liux1299/Software/snpEff/snpEff.jar"

VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.SNPs.private.vcf.gz"

# In dir: ~/Software/snpEff
# Build custom database since VCF is relative to parts reference
# Followed tutorial: https://pcingola.github.io/SnpEff/se_build_db_gff_gtf/
# Prepare directory
cd ~/Software/snpEff
mkdir data/Hordeum_vulgare_v3_parts

# Copy gff to directory
#cp /panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/gene_annotation/Hv_Morex.pgsb.Jul2020.sorted.parts.gff3 data/Hordeum_vulgare_v3_parts/genes.gff
# Add cds fasta sequence to end of copied gff file per SnpEff's documentation:
# https://pcingola.github.io/SnpEff/se_build_db_gff_gtf/#gff
# Add special comment line then fasta sequence
#echo "##FASTA" >> data/Hordeum_vulgare_v3_parts/genes.gff
#cat /panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/gene_annotation/Hv_Morex.pgsb.Jul2020.cds.fa >> data/Hordeum_vulgare_v3_parts/genes.gff

# Add entry to config file
echo "Hordeum_vulgare_v3_parts.genome: Hordeum_vulgare_v3_parts" >> snpEff.config
# Build database
java -jar ${SNPEFF_JAR} build -gtf22 Hordeum_vulgare_v3_parts

java -jar ${SNPEFF_JAR} Hordeum_vulgare ${VCF}
```
