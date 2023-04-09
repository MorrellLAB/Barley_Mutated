#!/usr/bin/env python
# Write by Li Lei 2021/06/15
# Modified by Chaochih Liu 2023/04/03 to handle slightly different file format
"""Decide which SNP is deleterious by BAD_Mutations. Since compile in BAD_Mutation is not correct if we
used 61 aligned genome, so we have to switch to the heuristics approach according to Kono et al., 2017 MBE
paper and Chun and Fay 2009 Genome Res paper.

Works with BAD_Mutations predictions compiled report combined with VeP report
(see script merge_bad_mut_and_vep.R to generate input file)
"""

import sys
import os

# User provided input arguments
# File should contain BAD_Mutations predictions merged with VeP report that contains amino acid info
# Assumes no "missing data" - some SNPs are not predicted by BAD_Mutations because of alignment problems
fp = os.path.expanduser(sys.argv[1])
# Test filepath
#fp = os.path.expanduser("~/Dropbox/Projects/Barley_Mutated/analyses/bad_mutations/predictions/mut_snps_private_Combined_Report_with_VeP.txt")
# Criteria
# Example: 10
min_seq = int(sys.argv[2])
# Example: 1
max_constraint = float(sys.argv[3])
# Example: 0.05
p_cutoff = float(sys.argv[4])
# We predicted on XX codons (nonsynonymous SNPs), so that becomes the basis of our multiple testing correction
# Number of codons (nonsynonymous SNPs) tested (given as input to BAD_Mutations predict)
# Here, it is the number of nonsynonymous SNPs that intersect with primary transcripts
# Example 611 nonsynonymous SNPs
num_codons = int(sys.argv[5])
lrt_sig = p_cutoff/num_codons

with open(fp, 'r') as f:
    for index, line in enumerate(f):
        # First line is header line
        if line.startswith('VariantID'):
            header = line.strip().split('\t')
            header.append('Del_vs_Tol')
            print('\t'.join(header))
        else:
            snp = line.strip().split('\t')
            # The relevant column indices after BAD_Mutations predictions are combined with VeP report:
            #   Alignment: 10
            #   MaskedConstraint: 12
            #   MaskedP.value: 13
            #   SeqCount: 9
            #   ReferenceAA: 11
            #   Amino_acids: 24
            aln = snp[10]
            masked_constraint = float(snp[12])
            masked_pval = float(snp[13])
            seq_count = int(snp[9])
            ref_aa_column = snp[11]
            amino_acids = snp[24].split('/')
            ref_aa = amino_acids[0]
            alt_aa = amino_acids[1]
            # Extra check to make sure ref amino acids match in merged input file
            if ref_aa_column == ref_aa:
                rn = aln.count(ref_aa)
                an = aln.count(alt_aa)
                # (min([rn,an]) == 0) -> either the alt or ref allele was not seen in any of the other species
                if (masked_pval < lrt_sig) and (seq_count >= min_seq) and (masked_constraint < max_constraint) and (min([rn,an]) == 0):
                    snp.append("Deleterious")
                else:
                    snp.append("Tolerated")
            else:
                raise Exception("Reference amino acids don't match for SNP: {}".format(snp[0]))
            print ('\t'.join(snp))
