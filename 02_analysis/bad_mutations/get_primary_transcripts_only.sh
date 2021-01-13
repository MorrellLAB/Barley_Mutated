#!/bin/bash

set -e
set -o pipefail

# After downloading the primary transcripts only files from Phytozome13, this script
#	generates a list of unique primary transcripts for bad_mutations predict

# User provided input arguments
# List of files that contain transcripts only
LIST_OF_FILES=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/phytozome13_download/phytozome/Hvulgare/r1/annotation/barley_primary_transcripts_only_files_list.txt
# Full filepath to output directory
OUT_DIR=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations
# Prefix of final output file
OUT_PREFIX=hvulgare

#------------------------------
# Prepare array
LIST_ARR=($(cat ${LIST_OF_FILES}))

# For each file, pull out primary transcript ID only, save to temp file
for i in ${LIST_ARR[@]}
do
	curr_prefix=$(basename ${i} .fa.gz)
	zgrep ">" ${i} | cut -d' ' -f 1 | sed -e 's,>,,' > ${OUT_DIR}/temp_${curr_prefix}.txt
done

# Combine and sort list of primary transcripts
cat ${OUT_DIR}/temp_*primaryTranscriptOnly.txt | sort -Vu > ${OUT_DIR}/${OUT_PREFIX}_primary_transcripts_only.txt

# Cleanup intermediate files
rm ${OUT_DIR}/temp_*primaryTranscriptOnly.txt
