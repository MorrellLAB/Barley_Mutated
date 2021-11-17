#!/bin/bash

# This script takes in a a Slurm job ID number and generates a list of
#   bad mutations predict re-run array indices
# Example: For a list of numbers (one per line), 1 2 3 5 6 7 9 10 11 12
#   this script collapses them into format 1-3,5-7,9,10-12

# Usage: ./get_re-run_array_indices.sh [xxxxxx]

# Where: [xxxxxx] is the slurm JobID number

set -e
set -o pipefail

# User provided input arguments
JOB_ID="$1"

#------------------------------
# Generate a list of timeout/failed bad_mutations predict array indices to re-run
# Grep invert match multiple strings: "JobID" "-----" "extern" "batch" "COMPLETED"
# Collapse consecutive numbers into ranges using awk
sacct -j ${JOB_ID} | grep -v 'JobID\|-----\|extern\|batch\|COMPLETED' | \
    awk '{ print $1 }' | sed -e "s/${JOB_ID}_//" | sort -V | \
    awk '
        function output() { print start (prev == start ? "" : "-"prev) }
        NR == 1 {start = prev = $1; next}
        $1 > prev+1 {output(); start = $1}
        {prev = $1}
        END {output()}
    ' | tr '\n' ',' | sed '$s/,$/\n/'