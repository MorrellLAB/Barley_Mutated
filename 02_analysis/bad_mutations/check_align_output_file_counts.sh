#!/bin/bash

set -e
set -o pipefail

# This script checks that we have the expected number of *.fa and *.tree files written to the MSA_Output directory.

# User provided input arguments
# We'll need the Fasta list of lists defined in the bad_mut_align.sh script
# Fasta list of lists (e.g., all_cds_hvulgare_list_of_lists.txt)
FASTA_LIST_OF_LISTS="$1"
# Path to MSA_output directory
MSA_OUTPUT_DIR="$2"

#------------------------------
# Make sure we are in the MSA_Output directory when we run this!
cd ${MSA_OUTPUT_DIR}

# Start from clean list since we will be appending to the output files
# If they already exist, delete
if [ -f ${MSA_OUTPUT_DIR}/temp_msa_output_problem_dirs.txt ]; then
    # File exists, remove before proceeding
    rm ${MSA_OUTPUT_DIR}/temp_msa_output_problem_dirs.txt
fi

if [ -f ${MSA_OUTPUT_DIR}/temp_msa_output_missing_transcripts_log_files.txt ]; then
    # File exists, remove before proceeding
    rm ${MSA_OUTPUT_DIR}/temp_msa_output_missing_transcripts_log_files.txt
fi

if [ -f ${MSA_OUTPUT_DIR}/temp_msa_output_evalue_error_log_files.txt ]; then
    rm ${MSA_OUTPUT_DIR}/temp_msa_output_evalue_error_log_files.txt
fi

if [ -f ${MSA_OUTPUT_DIR}/temp_msa_output_unexpected_symbol_error_log_files.txt ]; then
    rm ${MSA_OUTPUT_DIR}/temp_msa_output_unexpected_symbol_error_log_files.txt
fi

if [ -f ${MSA_OUTPUT_DIR}/temp_msa_output_other_error_log_files.txt ]; then
    rm ${MSA_OUTPUT_DIR}/temp_msa_output_other_error_log_files.txt
fi

# Figure out which ones need to be investigated/re-run
for i in $(cat ${FASTA_LIST_OF_LISTS})
do
    msa_subdir_name=$(basename ${i} .txt)
    fa_dir_name=$(dirname ${i})
    msa_output_dir=$(pwd)
    # The .tree files get written last and should be the same as the number of transcripts
    expected_tree_count=$(wc -l ${i} | cut -d' ' -f 1)
    actual_tree_count=$(ls ${msa_subdir_name}/*.tree | wc -l)
    if [ ${actual_tree_count} -ne ${expected_tree_count} ]; then
        echo "${msa_subdir_name}: ${actual_tree_count}"
        echo "${msa_subdir_name}: ${actual_tree_count}" >> ${MSA_OUTPUT_DIR}/temp_msa_output_problem_dirs.txt
        # Figure out which transcripts didn't get run to completion
        # Since we are appending, make sure we start from a clean list every time
        if [ -f ${msa_subdir_name}/all_log_files/temp_completed_tree_list.txt ]; then
            # File exists, remove before proceeding
            rm ${msa_subdir_name}/all_log_files/temp_completed_tree_list.txt
        fi
        for t in $(ls ${msa_subdir_name}/*.tree)
        do
            basename ${t} .tree >> ${msa_subdir_name}/all_log_files/temp_completed_tree_list.txt
        done
        # Again, start from a clean list in case we re-run this same command
        if [ -f ${msa_subdir_name}/all_log_files/temp_missing_transcripts.txt ]; then
            rm ${msa_subdir_name}/all_log_files/temp_missing_transcripts.txt
        fi
        # Create a list of log files for transcripts that didn't get run to completion
        grep -vf ${msa_subdir_name}/all_log_files/temp_completed_tree_list.txt ${i} >> ${msa_subdir_name}/all_log_files/temp_missing_transcripts.txt
        if [ -f ${msa_subdir_name}/all_log_files/temp_missing_transcripts_log_files.txt ]; then
            rm ${msa_subdir_name}/all_log_files/temp_missing_transcripts_log_files.txt
        fi
        for lf in $(cat ${msa_subdir_name}/all_log_files/temp_missing_transcripts.txt)
        do
            tname=$(basename ${lf} .fa)
            find ${msa_output_dir}/${msa_subdir_name}/all_log_files -name "${tname}*" >> ${msa_subdir_name}/all_log_files/temp_missing_transcripts_log_files.txt
        done
        echo "Printing list of filepaths to transcript log files and saving to file..."
        cat ${msa_subdir_name}/all_log_files/temp_missing_transcripts_log_files.txt >> ${MSA_OUTPUT_DIR}/temp_msa_output_missing_transcripts_log_files.txt
        echo "Missing transcripts are available here: ${msa_output_dir}/${msa_subdir_name}/all_log_files/temp_missing_transcripts.txt"
        printf "\n"
    fi
done

echo "List of all log files associated with missing transcripts are available here: ${MSA_OUTPUT_DIR}/temp_msa_output_missing_transcripts_log_files.txt"
# Pull only uniq lines
sort -uV ${MSA_OUTPUT_DIR}/temp_msa_output_missing_transcripts_log_files.txt > ${MSA_OUTPUT_DIR}/temp_msa_output_missing_transcripts_log_files_uniq.txt
mv ${MSA_OUTPUT_DIR}/temp_msa_output_missing_transcripts_log_files_uniq.txt ${MSA_OUTPUT_DIR}/temp_msa_output_missing_transcripts_log_files.txt

echo "Separating log files into categories: 1) e-value threshold related, 2) Unexpected symbol ('*',)"

for f in $(cat ${MSA_OUTPUT_DIR}/temp_msa_output_missing_transcripts_log_files.txt)
do
    if grep -wq "Could not find any BLAST hits! Try raising the E-value threshold for homology" ${f}
    then
        # Pattern if found
        echo ${f} >> ${MSA_OUTPUT_DIR}/temp_msa_output_evalue_error_log_files.txt
    elif grep -wq "Unexpected symbol" ${f}
    then
        # Pattern if found
        # Full pattern is "Unexpected symbol ('*',) in file of datatype PROTEIN"
        echo ${f} >> ${MSA_OUTPUT_DIR}/temp_msa_output_unexpected_symbol_error_log_files.txt
    else
        # Some other error that is not either of the above
        echo ${f} >> ${MSA_OUTPUT_DIR}/temp_msa_output_other_error_log_files.txt
    fi
done

# The above leads to duplicate files being appended to the final file, so we'll pull just uniq lines
# Evalue errors
sort -uV ${MSA_OUTPUT_DIR}/temp_msa_output_evalue_error_log_files.txt > ${MSA_OUTPUT_DIR}/temp_msa_output_evalue_error_log_files_uniq.txt
mv ${MSA_OUTPUT_DIR}/temp_msa_output_evalue_error_log_files_uniq.txt ${MSA_OUTPUT_DIR}/temp_msa_output_evalue_error_log_files.txt
# Unexpected symbol errors
sort -uV ${MSA_OUTPUT_DIR}/temp_msa_output_unexpected_symbol_error_log_files.txt > ${MSA_OUTPUT_DIR}/temp_msa_output_unexpected_symbol_error_log_files_uniq.txt
mv ${MSA_OUTPUT_DIR}/temp_msa_output_unexpected_symbol_error_log_files_uniq.txt ${MSA_OUTPUT_DIR}/temp_msa_output_unexpected_symbol_error_log_files.txt
# Other errors
sort -uV ${MSA_OUTPUT_DIR}/temp_msa_output_other_error_log_files.txt > ${MSA_OUTPUT_DIR}/temp_msa_output_other_error_log_files_uniq.txt
mv ${MSA_OUTPUT_DIR}/temp_msa_output_other_error_log_files_uniq.txt ${MSA_OUTPUT_DIR}/temp_msa_output_other_error_log_files.txt

echo "Lists of log files separated by error type are located in the following files:"
echo "${MSA_OUTPUT_DIR}/temp_msa_output_evalue_error_log_files.txt"
echo "${MSA_OUTPUT_DIR}/temp_msa_output_unexpected_symbol_error_log_files.txt"
echo "${MSA_OUTPUT_DIR}/temp_msa_output_other_error_log_files.txt"
