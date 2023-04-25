#!/bin/bash

# Log of commands run to check to figure out where errors are occurring

# User provided input arguments
# Full path to a list of lists (to utilize GNU parallel and job arrays)
# List of lists here should only include FASTA files that had an alignment
#	(i.e., .fa and .tree files were generated) in the previous align step
FASTA_LIST_OF_LISTS="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/final_lists/all_cds_hvulgare_list_of_lists.txt"

# Full path to a list of MSA_Output directories that contain *.fa and *.tree files
MSA_DIR_LIST="/panfs/jay/groups/9/morrellp/shared/Projects/WBDC_inversions/bad_mutations/results/MSA_output/MSA_out_dir_list_of_dir.txt"

#---------------------
for cds_list in $(cat ${FASTA_LIST_OF_LISTS})
do
    echo $cds_list
    for fasta_file in $(cat $cds_list)
    do
        echo $fasta_file
        # Example: hvulgare_cds_list-000
        curr_dir=$(basename $fasta_file .txt)
        # Example: /panfs/jay/groups/9/morrellp/shared/Projects/WBDC_inversions/bad_mutations/results/MSA_output/hvulgare_cds_list-000
        # Contains the .tree and .fasta files output from the align step
        curr_msa_dir=$(grep -w $curr_dir ${MSA_DIR_LIST})
        # Example: HORVU.MOREX.r3.1HG0000010.1
        curr_cds_name=$(basename $fasta_file .fasta)
        # Check if transcript is in the right list

    done
done
