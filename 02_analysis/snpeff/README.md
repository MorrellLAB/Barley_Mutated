# SnpEff

Run SnpEff as a check/comparison to premature stops identifed by VeP and Annovar.

Run SnpEff on mutated private SNPs.

Convert parts to pseudo positions.

```bash
# ~/Projects/Mutant_Barley/de_novo_vcfs
module load python3/3.8.3_anaconda2020.07_mamba

~/GitHub/File_Conversions/Barley_Parts_to_Pseudomolecules.py --vcf mut8_and_3mut10xGenomics.SNPs.private.vcf.gz morex_v3 > mut8_and_3mut10xGenomics.SNPs.private.pseudo_pos.vcf

~/GitHub/File_Conversions/Barley_Parts_to_Pseudomolecules.py --vcf mut8_and_3mut10xGenomics.INDELs.private.vcf.gz morex_v3 > mut8_and_3mut10xGenomics.INDELs.private.pseudo_pos.vcf
```

Use Morex v3 Hordeum_vulgare database from SnpEff.

```bash
# In dir: ~/Software/snpEff
# Dependencies
module load java/openjdk-17.0.2
SNPEFF_JAR="/panfs/jay/groups/9/morrellp/liux1299/Software/snpEff/snpEff.jar"

# SNPs
VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.SNPs.private.pseudo_pos.vcf"
PREFIX="mut8_and_3mut10xGenomics.SNPs.private.pseudo_pos"
# Indels
VCF_INDELs="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.INDELs.private.pseudo_pos.vcf"
PREFIX_INDELs="mut8_and_3mut10xGenomics.INDELs.private.pseudo_pos"
# Output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/snpeff"

# Annotate VCF
# SNPs
java -jar ${SNPEFF_JAR} -v -stats ${OUT_DIR}/${PREFIX}.html Hordeum_vulgare ${VCF} > ${OUT_DIR}/${PREFIX}.ann.vcf
# indels
java -jar ${SNPEFF_JAR} -v -stats ${OUT_DIR}/${PREFIX_INDELs}.html Hordeum_vulgare ${VCF_INDELs} > ${OUT_DIR}/${PREFIX_INDELs}.ann.vcf
```

Reformat output file for comparisons.

```bash
module load bcftools/1.10.2

# Test
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' mut8_and_3mut10xGenomics.SNPs.private.pseudo_pos.ann.vcf | tr '|' '\t'
```
