#!/usr/bin/env python3
"""Using the regions.bed.gz output from mosdepth, extract
high depth regions. This is defined as regions with coverage
greater than twice the standard deviation above the mean
coverage of the sample.

Usage: ./extract_high_depth_regions.py [regions.bed.gz] [sample_name] [out_dir]

Where:
1) [regions.bed.gz] is a file output from mosdepth for one sample
2) [sample_name] is the sample name used for the output file
3) [out_dir] is the path to the directory to store the output files
"""

import sys
import os
import pandas as pd
import numpy as np

# User provided input argumnets
# regions.bed.gz file
sample_bed = os.path.expanduser(sys.argv[1])
sample_name = sys.argv[2]
# Strip trailing slash so output path doesn't get messed up
out_dir = os.path.expanduser(sys.argv[3]).rstrip("/")

# Print usage message if no input arguments are provided
if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(0)

# Read mosdepth regions.bed.gz file
mosdepth_bed = pd.read_csv(
               sample_bed,
               compression = 'gzip',
               sep = '\t',
               names = ['chrom','start', 'end', 'depth'])

# Do some calculations
mean_depth = mosdepth_bed.depth.mean()
std_depth = mosdepth_bed.depth.std()
high_cov_bed = mosdepth_bed.loc[mosdepth_bed.depth > (mean_depth + 2*std_depth)]

out_file = os.path.expanduser(out_dir) + '/' + sample_name + '.high_cov.bed'
high_cov_bed.to_csv(out_file, sep = '\t', header = False, index = False)
