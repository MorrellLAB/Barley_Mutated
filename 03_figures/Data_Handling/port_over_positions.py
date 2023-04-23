#!/usr/bin/env python3

"""This script takes in 1) two VCF files (BOPA and 9k) containing the physical
positions of a set of SNPs and 2) a VCF file that needs the positions updated.
The output VCF file will contain updated physical positions. Known limitations:
The VCF should not be extremely large in file size.

This script outputs files in the same directory as the [VCF_to_update] file.

Usage: ./port_over_positions.py [BOPA_VCF] [9k_VCF] [VCF_to_update]
"""

import os
import sys
import gzip

if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(0)


def read_vcf(filename):
    """Read in a VCF file and store in dictionary."""
    header_lines = []
    vcf_dict = {}
    if '.gz' in filename:
        with gzip.open(filename, "rt") as file:
            for line in file:
                if line.startswith('#'):
                    header_lines.append(line.strip())
                else:
                    tmp = line.strip().split('\t')
                    vcf_dict[tmp[2]] = tmp
    else:
        with open(filename, "rt") as file:
            for line in file:
                if line.startswith('#'):
                    header_lines.append(line.strip())
                else:
                    tmp = line.strip().split('\t')
                    vcf_dict[tmp[2]] = tmp
    return header_lines, vcf_dict


def main(VCF_BOPA, VCF_9K, VCF_TO_UPDATE):
    """Driver function."""
    # Read in VCF files
    header_bopa, vcf_bopa = read_vcf(os.path.expanduser(VCF_BOPA))
    header_9k, vcf_9k = read_vcf(os.path.expanduser(VCF_9K))
    header_update, vcf_update = read_vcf(os.path.expanduser(VCF_TO_UPDATE))
    updated_dict = {}
    no_match = {}
    # Update chr and position
    for key in vcf_update:
        if key in vcf_bopa and key in vcf_9k:
            # Use position in 9k set
            updated_dict[key] = vcf_9k[key][0:2] + vcf_update[key][2:]
        elif key in vcf_bopa and key not in vcf_9k:
            # SNP only exists in BOPA set
            updated_dict[key] = vcf_bopa[key][0:2] + vcf_update[key][2:]
        elif key in vcf_9k and key not in vcf_bopa:
            # SNP only exists in 9k set
            updated_dict[key] = vcf_9k[key][0:2] + vcf_update[key][2:]
        else:
            # SNP is not in either set
            no_match[key] = vcf_update[key]
    # Generate output file prefix
    out_dir = os.path.dirname(VCF_TO_UPDATE)
    out_prefix = os.path.basename(os.path.splitext(VCF_TO_UPDATE)[0])
    out_file = out_dir + '/' + out_prefix + '_physPos.vcf'
    no_match_out = out_dir + '/' + out_prefix + '_noMatch.vcf'
    # Save SNPs to file
    if os.path.isfile(out_file):
        # Start from a clean file
        os.remove(out_file)
        os.remove(no_match_out)
        # Save header lines to output files
        for i in header_update:
            with open(out_file, 'a') as out:
                out.write(i + '\n')
        for i in header_update:
            with open(no_match_out, 'a') as out:
                out.write(i + '\n')
        # Save current output to file
        for key in updated_dict:
            with open(out_file, 'a') as out:
                out.write('\t'.join(updated_dict[key]) + '\n')
        # Save SNPs not in BOPA or 9k set
        for key in no_match:
            with open(no_match_out, 'a') as out:
                out.write('\t'.join(no_match[key]) + '\n')
    else:
        # Save header lines to output files
        for i in header_update:
            with open(out_file, 'a') as out:
                out.write(i + '\n')
        for i in header_update:
            with open(no_match_out, 'a') as out:
                out.write(i + '\n')
        # Save current output to file
        for key in updated_dict:
            with open(out_file, 'a') as out:
                out.write('\t'.join(updated_dict[key]) + '\n')
        # Save SNPs not in BOPA or 9k set
        for key in no_match:
            with open(no_match_out, 'a') as out:
                out.write('\t'.join(no_match[key]) + '\n')


main(sys.argv[1], sys.argv[2], sys.argv[3]) # Run the program
