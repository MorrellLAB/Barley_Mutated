#!/usr/bin/env python3
"""Script that splits chromosomes listed in .fasta.fai files into N parts and
formats the output into GATK's .intervals file format: <chr>:<start>-<stop>
Only the <chr> part is strictly required, see here for more details on
format:
https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists

Currently, chrUn does not get split and the script will need to be modified if there are scaffolds present.

Usage: ./split_chrom_fai.py [fai_filepath] [num_parts] > /path/to/chromosomes_split.intervals

Where:

1) [fai_filepath] is the path to the .fasta.fai file
2) [num_parts] is the number of parts we want to split each
chromosome in the .fasta.fai file into (integer).
"""

import sys
import os
import math

# Print usage message if no input arguments are provided
if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(0)

def read_fai(fai_filepath):
    """Read in .fasta.fai file and store in dictionary."""
    fai_dict = {}
    with open(fai_filepath, 'rt') as f:
        for l in f:
            lrecord = l.strip().split()
            chr_name = lrecord[0]
            fai_dict[chr_name] = lrecord
    return fai_dict


def main(fai_filepath, num_parts):
    """Driver function."""
    fai_dict = read_fai(os.path.expanduser(fai_filepath))
    # Split chromosomes into N parts
    split_intervals = []
    for key in fai_dict:
        if key != "chrUn":
            curr_chr = fai_dict[key][0]
            chr_length = fai_dict[key][1]
            # Divide length of chr into n parts
            part_size = math.ceil(int(chr_length) / int(num_parts))
            # Note: second number is exclusive so add 1
            nparts = list(range(1, int(num_parts) + 1))
            # Generate new chr intervals
            start_pos = "1"
            end_pos = chr_length
            for p in nparts:
                if p != nparts[-1]:
                    # If we are not on the last loop, generate current end pos
                    curr_end_pos = int(start_pos) + (int(part_size) - 1)
                    new_interval = curr_chr + ":" + str(start_pos) + "-" + str(curr_end_pos)
                    split_intervals.append(new_interval)
                    # Update start position for next iteration
                    start_pos = curr_end_pos + 1
                elif p == nparts[-1]:
                    # If we are on the last loop, use the end_pos
                    new_interval = curr_chr + ":" + str(start_pos) + "-" + str(end_pos)
                    split_intervals.append(new_interval)
        elif key == "chrUn":
            # Don't need to split
            split_intervals.append(fai_dict[key][0])

    # Print new split intervals to stdout
    print('\n'.join(split_intervals))


main(sys.argv[1], sys.argv[2])
