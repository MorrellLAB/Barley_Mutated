#!/bin/bash

set -e
set -o pipefail

# This script prepares list of lists of fasta files for a seperate run where we try increasing the e-values for trasncripts that returned this error.

# User provided input arguments
# List of log files for transcripts that have the e-value error
evalue_err_list="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/MSA_output/temp_msa_output_evalue_error_log_files.txt"
# Filepath to our original align lists directory
align_lists_dir="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/align_lists"
# Filepath to our output directory that'll store the new align lists
out_dir_align_lists="/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/align_lists_increased_e-val"

#----------------------
# Make output directory
mkdir -p ${out_dir_align_lists}
# Since we are apending files, start from clean slate
# Delete output directory contents before starting each time
if [ "$(ls -A ${out_dir_align_lists})" ]
then
    # Directory is not empty
    echo "Output directory has files, clearing to start from a clean slate since we are appending to files."
    rm ${out_dir_align_lists}/*.txt
fi

# Store list of unique directories so our list of lists gets structured
#   in the same way as the original run of bad_mut_align.sh
uniq_dir_arr=($(sed -e 's,\(.*\)/,\1\t,' ${evalue_err_list} | cut -f 1 | sort -uV | sed -e 's,/all_log_files,,' | sed -e 's,\(.*\)/,\1\t,' | cut -f 2))

for d in ${uniq_dir_arr[@]}
do
    echo ${d}
    for i in $(grep ${d} ${evalue_err_list})
    do
        echo ${i}
        curr_transcript=$(basename ${i} .log)
        # Pull out fasta file path
        grep -w "${curr_transcript}" ${align_lists_dir}/hvulgare_cds_list-*.txt | cut -d":" -f 2 >> ${out_dir_align_lists}/${d}.txt
    done
done
