#!/usr/bin/env python3
"""Convert from ANNOVAR output to a unified table, in the same form as the
output from SNP_Effect_Predictor.py:
    SNP ID
    Chromosome
    Position (1-based)
    Silent or not?
    Transcript ID
    Codon position
    Reference allele
    Alternate allele
    AA1
    AA2
    Residue number
"""

# This script originated from Tom Kono and was modified by Li Lei
#   CL updated it to work with formats output from ANNOVAR's latest version of their tools installed as of 2021-11-18
#   Tom's version: https://github.com/MorrellLAB/File_Conversions/blob/master/ANNOVAR_To_Effects.py
#   Li's version: https://github.com/lilei1/MBE_samples/blob/master/Script/ANNOVA_To_effects.py

import sys
import os

all_effects = os.path.expanduser(sys.argv[1])
exon_effects = os.path.expanduser(sys.argv[2])

coding_variants = {}
with open(exon_effects, 'r') as f:
    for line in f:
        tmp = line.strip().split('\t')
        if tmp[1] == 'synonymous SNV':
            silent = 'Yes'
        else:
            silent = 'No'
        ann = tmp[2].split(':')
        txid = ann[1]
        aa1 = ann[4][-2]
        aa2 = ann[4][2]
        if tmp[1] == 'stopgain':
            aa2 = '*'
        if tmp[1] == 'stoploss':
            aa1 = '*'
        aapos = ann[4][3:-2]
        pos = tmp[4]
        chrom = tmp[3]
        snpid = '.'.join([chrom,pos])
        ref_allele = tmp[6]
        alt_allele = tmp[7]
        cdspos = ann[3][3:-1]
        codonpos = str(3 - (int(cdspos) % 3))
        coding_variants[snpid] = (
            chrom,
            pos,
            silent,
            txid,
            codonpos,
            ref_allele,
            alt_allele,
            aa1,
            aa2,
            aapos)

print ('\t'.join(
    [
        'SNP_ID',
        'Chromosome',
        'Position',
        'Silent',
        'Transcript_ID',
        'Codon_Position',
        'Ref_Base',
        'Alt_Base',
        'AA1',
        'AA2',
        'CDS_Pos'
    ]))

with open(all_effects, 'r') as f:
    for line in f:
        tmp = line.strip().split('\t')
        chrom = tmp[2]
        pos = tmp[3]
        snpid = '.'.join([chrom,pos])

        if tmp[0] == 'exonic':
            chrom, pos, silent, txid, codonpos, ref_allele, alt_allele, aa1, aa2, aapos = coding_variants[snpid]
        else:
            silent = 'Yes'
            txid = '-'
            codonpos = '-'
            ref_allele = tmp[5]
            alt_allele = tmp[6]
            aa1 = '-'
            aa2 = '-'
            aapos = '-'
        print ('\t'.join(
            [
                snpid,
                chrom,
                pos,
                silent,
                txid,
                codonpos,
                ref_allele,
                alt_allele,
                aa1,
                aa2,
                aapos
            ]))
