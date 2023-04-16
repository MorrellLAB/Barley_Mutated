#!/usr/bin/env python3
"""Convert the table with Deleterious vs Tolerated info from
barley parts positions to pseudomolecules positions. This makes
visualization more convenient. Two new columns "Chr_pseudo" and "Pos_pseudo"
will be appended to the end of each line with the converted chr name and positions.

Usage: ./del_vs_tol_barley_parts_to_pseudomolecules.py [del_vs_tol.txt] > output_del_vs_tol.with_pseudo_pos.txt

Where:
1. [del_vs_tol.txt] is a file containing "Deleterious" vs "Tolerated" annotations and barley parts positions
"""

import sys
import os

if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(0)

#   Store the lengths of the barley pseudomolecule parts in a dictionary. We only
#   really need to store the first parts, since we will just add the part2
#   positions to it.
# Define the length of the part1 pieces as a dictionary constant.
# Morex v3 sizes
PARTS_SIZES_V3 = {
    'chr1H': 206486643,
    'chr2H': 301293086,
    'chr3H': 267852507,
    'chr4H': 276149121,
    'chr5H': 204878572,
    'chr6H': 256319444,
    'chr7H': 328847192,
    'chrUn': 29110253,
    'EF115541.1': 136462,
    'AP017301.1': 525599
}

# User provided input arguments
# Input snp table example: mut_snps_private_deleterious_vs_tolerated.txt
snp_table_fp = os.path.expanduser(sys.argv[1])

with open(snp_table_fp, 'rt') as f:
    for line in f:
        if line.startswith("VariantID"):
            header = line.strip().split('\t')
            header.extend(["Chr_pseudo", "Pos_pseudo"])
            print('\t'.join(header))
        else:
            tmp = line.strip().split('\t')
            # Example column 16: chr1H_part1:102622702
            parts_chr_pos = tmp[16].split(':')
            # Modify the chromosome to not have the part
            chrom = parts_chr_pos[0].split('_')[0]
            # Convert parts to pseudo positions
            if '_part1' in parts_chr_pos[0]:
                # Position can remain as is, modify chr name to not have the part
                pseudo_chr_pos = [chrom, parts_chr_pos[1]]
            elif '_part2' in parts_chr_pos[0]:
                # If chr name has '_part2' in it, then we have to modify the chr name and the position
                offset = PARTS_SIZES_V3[chrom]
                newpos = str(int(parts_chr_pos[1]) + offset)
                # Use new position
                pseudo_chr_pos = [chrom, newpos]
            # Append modified values to the end of list
            tmp.extend(pseudo_chr_pos)
            print('\t'.join(tmp))
