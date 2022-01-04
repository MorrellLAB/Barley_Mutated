#!/bin/bash

# This script keeps a record of commands run to split VCF by sample

# Dependencies
module load bcftools/1.9
module load python3/3.8.3_anaconda2020.07_mamba

# Running commands from the following directory
cd /panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated/Filtered

# Heterozygous and Homozygous singletons VCF
# Pull out deletions only for verification
bcftools filter -i 'INFO/SVTYPE=="DEL"' mut_3_lines_filtered_singletons_only_annotated_DEL_de_novo_sites.vcf > deletions/mut_3_lines_filtered_singletons_only-all_DEL.vcf
# Pull out deletions that were scored
cd deletions
bcftools filter -i 'INFO/SVSupScore >= 0' mut_3_lines_filtered_singletons_only-all_DEL.vcf > mut_3_lines_filtered_singletons_only-scored_DEL.vcf
bcftools filter -e 'INFO/SVSupScore < 75.0' mut_3_lines_filtered_singletons_only-scored_DEL.vcf > mut_3_lines_filtered_singletons_only-scored_DEL_gte75Sup.vcf
bcftools filter -e 'INFO/SVSupScore < 95.0' mut_3_lines_filtered_singletons_only-scored_DEL.vcf > mut_3_lines_filtered_singletons_only-scored_DEL_gte95Sup.vcf
# Split by sample
bcftools view -s M01-3-3 mut_3_lines_filtered_singletons_only-scored_DEL_gte75Sup.vcf | bcftools filter -e 'GT="mis"' > M01_singletons_only-scored_DEL_gte75Sup.vcf
bcftools view -s M20-2-2 mut_3_lines_filtered_singletons_only-scored_DEL_gte75Sup.vcf | bcftools filter -e 'GT="mis"' > M20_singletons_only-scored_DEL_gte75Sup.vcf
bcftools view -s M29-2-2 mut_3_lines_filtered_singletons_only-scored_DEL_gte75Sup.vcf | bcftools filter -e 'GT="mis"' > M29_singletons_only-scored_DEL_gte75Sup.vcf
# Convert to BED format
# We found 10x Genomics lists end positions of SVs greater than 1 bp in the INFO field under END=. So, we'll need to use a custom script to do the VCF to BED file conversion for the deletions.
~/GitHub/Barley_Mutated/00_sequence_processing/01_snp_filtering/vcf_10x_genomics_to_bed.py M01_singletons_only-scored_DEL_gte75Sup.vcf > M01_singletons_only-scored_DEL_gte75Sup.bed
~/GitHub/Barley_Mutated/00_sequence_processing/01_snp_filtering/vcf_10x_genomics_to_bed.py M20_singletons_only-scored_DEL_gte75Sup.vcf > M20_singletons_only-scored_DEL_gte75Sup.bed
~/GitHub/Barley_Mutated/00_sequence_processing/01_snp_filtering/vcf_10x_genomics_to_bed.py M29_singletons_only-scored_DEL_gte75Sup.vcf > M29_singletons_only-scored_DEL_gte75Sup.bed

# Homozygous singletons VCF
cd /panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated/Filtered
# Pull out homozygous singletons
bcftools view -i 'GT[*]="hom"' mut_3_lines_filtered_singletons_only_annotated_DEL_de_novo_sites.vcf > mut_3_lines_filtered_hom_singletons_only_annotated_DEL_de_novo_sites.vcf
# Split by sample for counting
# M01
bcftools filter -e 'INFO/SVSupScore < 75.0' mut_3_lines_filtered_hom_singletons_only_annotated_DEL_de_novo_sites.vcf | bcftools view -s M01-3-3 | bcftools filter -e 'GT="mis"' > by_sample_svplaudit_scored/M01-3-3_hom_singletons_only_annotated_DEL_gte_75.vcf
bcftools filter -e 'INFO/SVSupScore < 75.0' mut_3_lines_filtered_singletons_only_annotated_DEL_de_novo_sites.vcf | bcftools view -s M01-3-3 | bcftools filter -e 'GT="mis"' > by_sample_svplaudit_scored/M01-3-3_singletons_only_annotated_DEL_gte_75.vcf
# M20
bcftools filter -e 'INFO/SVSupScore < 75.0' mut_3_lines_filtered_hom_singletons_only_annotated_DEL_de_novo_sites.vcf | bcftools view -s M20-2-2 | bcftools filter -e 'GT="mis"' > by_sample_svplaudit_scored/M20-2-2_hom_singletons_only_annotated_DEL_gte_75.vcf
bcftools filter -e 'INFO/SVSupScore < 75.0' mut_3_lines_filtered_singletons_only_annotated_DEL_de_novo_sites.vcf | bcftools view -s M20-2-2 | bcftools filter -e 'GT="mis"' > by_sample_svplaudit_scored/M20-2-2_singletons_only_annotated_DEL_gte_75.vcf
# M29
bcftools filter -e 'INFO/SVSupScore < 75.0' mut_3_lines_filtered_hom_singletons_only_annotated_DEL_de_novo_sites.vcf | bcftools view -s M29-2-2 | bcftools filter -e 'GT="mis"' > by_sample_svplaudit_scored/M29-2-2_hom_singletons_only_annotated_DEL_gte_75.vcf
bcftools filter -e 'INFO/SVSupScore < 75.0' mut_3_lines_filtered_singletons_only_annotated_DEL_de_novo_sites.vcf | bcftools view -s M29-2-2 | bcftools filter -e 'GT="mis"' > by_sample_svplaudit_scored/M29-2-2_singletons_only_annotated_DEL_gte_75.vcf

# Count homozygous singletons
cd by_sample_svplaudit_scored
# M01
~/GitHub/Barley_Mutated/01_snp_filtering/count_SV_events.sh M01-3-3_singletons_only_annotated_DEL_gte_75.vcf
~/GitHub/Barley_Mutated/01_snp_filtering/count_SV_events.sh M01-3-3_hom_singletons_only_annotated_DEL_gte_75.vcf
# M20
~/GitHub/Barley_Mutated/01_snp_filtering/count_SV_events.sh M20-2-2_singletons_only_annotated_DEL_gte_75.vcf
~/GitHub/Barley_Mutated/01_snp_filtering/count_SV_events.sh M20-2-2_hom_singletons_only_annotated_DEL_gte_75.vcf
# M29
~/GitHub/Barley_Mutated/01_snp_filtering/count_SV_events.sh M29-2-2_singletons_only_annotated_DEL_gte_75.vcf
~/GitHub/Barley_Mutated/01_snp_filtering/count_SV_events.sh M29-2-2_hom_singletons_only_annotated_DEL_gte_75.vcf
