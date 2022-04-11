#  Samplot on mutated lines of barley

This file reports steps and scripts to run [Samplot](https://github.com/ryanlayer/samplot) and upload data to Amazon Web Service (AWS) account to perform the scoring of the images using [SV-Plaudit](https://github.com/jbelyeu/SV-plaudit) .

We dicided to focus our attention on SVs closed to genes. 

__To have a .bed file with gene interval positions increased by 20%__
```
# To select only the first 8 columns of gff3 file
cut -f1,2,3,4,5,6,7,8 Barley_Morex_V2_gene_annotation_PGSB.all.parts.gff3 > reduced.gff3

# To increase the intervals
bedtools slop -i reduced.gff3 -g barley.genome -b 0.2 -pct -header > add_ranges.bed
```
NOTE: barley.genome is a txt file in which the length of each chromosome or contig is reported to prevent the extension of intervals beyond chromosome boundaries [Bedtools-slop](https://bedtools.readthedocs.io/en/latest/content/tools/slop.html)

__To select only the SNPs inside the intervals__
```
bedtools intersect -a mut_3_lines_filtered_singletons_only_copy2.vcf -b add_ranges.bed -header > new_only_gene.vcf
```

__To run Samplot__

I installed Samplot on the ML environment on MSI at UMN.

NOTE: to call Samplot you have to check which anaconda environmment you are using.
```
# To check the anaconda environment
conda info --envs

# To change environment
source activate /home/morrellp/public/Software/anaconda3
```
We are lookig to __DELATIONS__.

```
# To run Samplot
samplot vcf --filter "SVTYPE == 'DEL'" --vcf /panfs/roc/scratch/gfrascar/barley_sodium/new_only_gene.vcf -n M01-3-3 M20-2-2 M29-2-2 Morex -d /panfs/roc/scratch/gfrascar/barley_sodium/png -O png --ped /panfs/roc/groups/9/morrellp/gfrascar/barley_sodium_png/mutant_ped.ped --dn_only -b /panfs/roc/groups/9/morrellp/shared/IGV/M01-3-3_phased_possorted_bam.bam /panfs/roc/groups/9/morrellp/shared/IGV/M20-2-2_phased_possorted_bam.bam /panfs/roc/groups/9/morrellp/shared/IGV/M29-2-2_phased_possorted_bam.bam /panfs/roc/groups/9/morrellp/shared/IGV/morex-sample2_phased_possorted_bam.bam -a > samplot_commands.sh
```
__IMPORTANT__
To run SV-Plaudit after Samplot the -a flag is mandatory. This option will generate .json.

__To run SV-Plaudit__

For creating the project in the AWS cloud you have to follow the guide in SV-Plaudit [GitHub page](https://github.com/jbelyeu/SV-plaudit)

__Annotations of the scores__

I used [VAtools](https://github.com/griffithlab/VAtools/issues/43)

REQUIRED:
1. A .tsv table in which are reported the __chromoosome__, __position__, and __value__ to be annotated.
2. A vcf file to be annotated

```
vcf-info-annotator -d "SV score in % for the Support value" -f Float -o mut_3_lines_filtered_singletons_only_annotated_DEL.vcf mut_3_lines_filtered_singletons_only.vcf test.tsv SVSupScore 
```
I annotated the value of the SUPPORT from the Summary report (that I preveusly downloade) 

The annotated vcf is:
```
/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated/Filtered/mut_3_lines_filtered_singletons_only_annotated_DEL.vcf
```




