#!/bin/bash

# This script counts the number of unique events per individual in each variant category
#   from 10x Genomics VCFs.

# Dependencies
module load bcftools/1.9

# User provided input arguments
# VCF filepath for a SINGLE sample
VCF=$(realpath "$1")

#----------------------
# SNPs
snps=$(bcftools filter -i 'TYPE=="snp"' "${VCF}" | grep -v "#" | wc -l)

# Insertion
# bcftools filter for some reason can't pull out TYPE=ins, so we'll use grep
ins=$(grep "TYPE=ins" "${VCF}" | grep -v "#" | wc -l)

# DELs
# DEL events that are NOT BND
dels_not_BND=$(bcftools filter -i 'INFO/SVTYPE=="DEL" | INFO/SVTYPE2=="DEL"' "${VCF}" | bcftools filter -e 'INFO/SVTYPE=="BND" | INFO/SVTYPE2=="BND"' | grep -v "#" | wc -l)
# DEL events that are also BND (need to count unique events here)
dels_BND=$(bcftools filter -i 'INFO/SVTYPE=="DEL" | INFO/SVTYPE2=="DEL"' "${VCF}" | bcftools filter -i 'INFO/SVTYPE=="BND" | INFO/SVTYPE2=="BND"' | grep -v "#" | cut -f 3 | sed -e 's,.$,,' | sort -uV | wc -l)
# Sum unique events together
total_dels=$(expr ${dels_not_BND} + ${dels_BND})

# DUP
# DUP events are separated into multiple VCF record entries, so the counting needs to account for this
# Remove last digit of calls and sort to get unique events
dup=$(bcftools filter -i 'INFO/SVTYPE=="DUP" | INFO/SVTYPE2=="DUP"' "${VCF}" | grep -v "#" | cut -f 3 | sed -e 's,.$,,' | sort -uV | wc -l)

# INV
# INV events are separated into multiple VCF record entries, so the counting needs to account for this
# Remove last digit of calls and sort to get unique events
#    Example: call_600_1 and call_600_2 are separate entries for the same INV event
inv=$(bcftools filter -i 'INFO/SVTYPE=="INV" | INFO/SVTYPE2=="INV"' "${VCF}" | grep -v "#" | cut -f 3 | sed -e 's,.$,,' | sort -uV | wc -l)
if [[ $? != 0 ]]; then
    # Command failed, returned nothing, set count to 0
    inv="0"
fi

# Total unique events in VCF
total_combined=$(expr ${snps} + ${ins} + ${total_dels} + ${dup} + ${inv})

# Print summary to stdout
echo "Number of SNPs: ${snps}"
echo "Number of insertions: ${ins}"
echo "Number of deletions: ${total_dels}"
echo "Number of duplications: ${dup}"
echo "Number of inversions: ${inv}"
echo "Total number of unique events: ${total_combined}"
