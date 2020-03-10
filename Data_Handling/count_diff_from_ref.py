#!/usr/bin/env python3

"""This script counts the number of differences from Morex given
T3 genotype data stored in a VCF format. Note: T3 genotype data was
downloaded in HapMap format and converted to VCF so the ref allele is
the major allele.

Usage: ./count_diff_from_ref.py [ref_name] [vcf_filename] [centromere_filename]
"""

import os
import sys
import pandas as pd

if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(0)


def read_vcf(filename):
    """Read in a VCF file and store in dictionary."""
    header_lines = []
    vcf_dict = {}
    with open(filename, "rt") as file:
        for line in file:
            if line.startswith('##'):
                header_lines.append(line.strip())
            elif line.startswith('#CHROM'):
                # Make program case insensitive, change to all uppercase
                tmp = line.strip().split('\t')
                header = []
                for i in tmp:
                    header.append(i.upper())
                vcf_dict[tmp[2]] = header
            else:
                tmp = line.strip().split('\t')
                vcf_dict[tmp[2]] = tmp
    return header_lines, vcf_dict


def read_centromere(filename):
    """Read in a file containing centromere positions for each
    chromosome."""
    cent_dict = {}
    with open(filename, "rt") as file:
        for line in file:
            # Skip header line
            if line.startswith('agp'):
                continue
            else:
                tmp = line.strip().split('\t')
                cent_dict[tmp[0]] = tmp[1]
    return cent_dict


def count_diffs(ref_name, ref_dat, line_name, line_dat, centromere):
    """Count the number of differences between the reference
    line and the current line we are processing. If either of
    the lines have missing data, do not count the line."""
    out_df = [ref_name, line_name, str(len(ref_dat))]
    # Total differences
    num_diff = 0
    num_missing = 0
    for i, elem in enumerate(ref_dat):
        # Do not count missing genotype
        if ref_dat[i] == "./." or line_dat[line_name][i] == "./.":
            num_missing += 1
        elif ref_dat[i] != line_dat[line_name][i]:
            num_diff += 1
            #print(ref_dat[i], line_dat[line_name][i])
    # Append total differences and n missing to list
    out_df.append(str(num_diff))
    out_df.append(str(num_missing))
    # Differences split by chromosome arms (L/R of centromere)
    for key in centromere:
        if key in list(set(line_dat['#CHROM'])):
            num_diff_left_arm = 0
            num_diff_right_arm = 0
            for i, elem in enumerate(ref_dat):
                if line_dat['#CHROM'][i] == key and int(line_dat['POS'][i]) <= int(centromere[key]):
                    if ref_dat[i] == "./." or line_dat[line_name][i] == "./.":
                        continue
                    elif ref_dat[i] != line_dat[line_name][i]:
                        num_diff_left_arm += 1
                elif line_dat['#CHROM'][i] == key and int(line_dat['POS'][i]) > int(centromere[key]):
                    if ref_dat[i] == "./." or line_dat[line_name][i] == "./.":
                        continue
                    elif ref_dat[i] != line_dat[line_name][i]:
                        num_diff_right_arm += 1
            # For each chr, append
            out_df.append(str(num_diff_left_arm))
            out_df.append(str(num_diff_right_arm))
    return out_df


def main(ref_line_name, vcf_filename, centromere_filename):
    """Driver function."""
    # Read in VCF file
    header, vcf = read_vcf(os.path.expanduser(vcf_filename))
    # Read in file containing centromere positions
    centromere_pos = read_centromere(os.path.expanduser(centromere_filename))
    # Store dict in pandas data frame
    pd_df = pd.DataFrame.from_dict(vcf, orient='index')
    # Make the first row the header
    new_header = pd_df.iloc[0]
    pd_df = pd_df[1:]
    pd_df.columns = new_header
    # Prepare list of line names
    line_names = []
    for i in new_header[9:]:
        if i == ref_line_name.upper():
            continue
        else:
            line_names.append(i)
    # Pull out reference line only
    ref_df = pd_df[ref_line_name.upper()]
    # Prepare output file header
    out_header = ['ref_line', 'line', 'total_markers', 'n_diffs', 'n_missing']
    cent_header = []
    # Prepare centromere header
    uniq_chrom = list(set(pd_df['#CHROM']))
    for key in centromere_pos:
        if key in uniq_chrom:
            cent_header.append(key + '_L')
            cent_header.append(key + '_R')
    for i in cent_header:
        out_header.append(i)
    # Count differences between lines
    out_diffs = []
    for i in line_names:
        cdf = count_diffs(ref_name=ref_line_name.upper(), ref_dat=ref_df, line_name=i, line_dat=pd_df, centromere=centromere_pos)
        # Save counts for each comparison
        out_diffs.append(cdf)
    # Print output file header line
    print('\t'.join(out_header))
    # Print counts to stdout
    for i, elem in enumerate(out_diffs):
        print('\t'.join(elem))


main(sys.argv[1], sys.argv[2], sys.argv[3]) # Run the program
