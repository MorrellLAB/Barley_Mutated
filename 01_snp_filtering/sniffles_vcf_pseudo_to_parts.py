#!/usr/bin/env python3
"""Convert from IPK full pseudomolecules to parts for Sniffles2 vcf formatted
SV calls (end position is stored in the INFO/END field and will also have its
position updated to the parts positions).

Important: currently only works for Morex v3.

Usage: ./sniffles_vcf_pseudo_to_parts.py [vcf_fp] > vcf_fp_parts.vcf

Where: [vcf_fp] is the Sniffles2 VCF file."""

import sys
import os

# Morex v3 sizes
parts_sizes = {
    'chr1H': 206486643,
    'chr2H': 301293086,
    'chr3H': 267852507,
    'chr4H': 276149121,
    'chr5H': 204878572,
    'chr6H': 256319444,
    'chr7H': 328847192,
    'chrUn': 29110253,
    'EF115541.1': 136462,
    'AP017301.1': 525599,
    'Pt': 999999999999
}

if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(0)

# User provided input arguments
vcf_fp = os.path.expanduser(sys.argv[1])

# If lines start with '#', print them
with open(vcf_fp, 'r') as f:
    for line in f:
        if line.startswith('#'):
            print(line.strip())
        else:
            tmp = line.strip().split()
            chrom = tmp[0]
            pos = int(tmp[1])
            # Check that the chromosomes are named as we are expecting
            if chrom not in parts_sizes:
                sys.stderr.write(chrom + ' not recognized. The chromosomes must be named like \'chr1H.\'\n')
            # Check the parts lengths. If the position is greater than the
            # part1 length, then subtract the part1 length, and change the
            # name to part2
            # Also passing 'chrUn' unchanged, as there is only 1 part
            limit = parts_sizes[chrom]
            if chrom == 'chrUn' or chrom == 'Pt' or chrom == 'EF115541.1' or chrom == 'AP017301.1':
                newchrom = chrom
                newpos = str(pos)
            elif pos > limit:
                newchrom = chrom + '_part2'
                newpos = str(pos - limit)
                # Check if INFO field contains the END= annotation
                #   If so, the END= position will need to be updated too
                #   Not all variants have this, it depends on the type of variant called
                if 'END=' in tmp[7]:
                    tmp_info = tmp[7].split(';')
                    # Identify index of END= field
                    end_idx = [i for i, s in enumerate(tmp_info) if 'END=' in s][0]
                    endpos = int(tmp_info[end_idx].split('=')[1])
                    newendpos = str(endpos - limit)
                    # Update END= field with new end position
                    tmp_info[end_idx] = '='.join(['END', newendpos])
                    # Rejoin INFO field
                    tmp[7] = ';'.join(tmp_info)
            else:
                newchrom = chrom + '_part1'
                newpos = str(pos)
            print('\t'.join([newchrom, newpos] + tmp[2:]))
