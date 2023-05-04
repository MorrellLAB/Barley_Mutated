8,10,14,
17,19,32-35,37-38,42-44,49,63,65-68,70-71,73-74,84,89,92-97,99,103-104,108,111,113-118,120-123,125,129-131,133-134,141-158,164,166,170-172,174-186,189-190,192-193

CURR_MSA_DIR=$(echo ${MSA_DIR_ARR[33]})

bad_mut_script="${BAD_MUT_SCRIPT}"
config_file="${CONFIG_FILE}"
fasta_list="${CURR_FASTA_LIST}"
msa_fasta="${MSA_FASTA_ARR[1]}"
curr_msa_dir="${CURR_MSA_DIR}"
subs_list="${OUT_DIR}/primary_transcript_subs_list.txt"
out_dir="${OUT_DIR}"




msa_name=$(basename ${msa_fasta} .fasta)
# Get msa_tree file
msa_tree="${curr_msa_dir}/${msa_name}.tree"
tree_name=$(basename ${msa_tree} .tree)

# Check that msa_fasta and msa_tree names match
if [[ ${msa_name} != ${tree_name} ]]
then
    echo "There is a mismatch between the MSA fasta ${msa_name} and tree file ${tree_name}, exiting..."
    exit 1
fi

# Get current list prefix
curr_list_prefix=$(basename ${fasta_list} .txt)
# Get fasta file that corresponds with MSA fasta file
fasta_file=$(grep -w ${msa_name} ${fasta_list})

# Get subs filename
subs_file=$(grep -w ${msa_name} ${subs_list})
subdir_name="${curr_list_prefix}"
# Check if out subdirectory exists, if not make it
mkdir -p ${out_dir}/${subdir_name} ${out_dir}/all_predict_log/${subdir_name} ${out_dir}/${subdir_name}/problematic_predictions
echo "Config file: ${config_file}"
echo "Fasta file: ${fasta_file}"
echo "MSA fasta: ${msa_fasta}"
echo "MSA tree: ${msa_tree}"
echo "Subs file: ${subs_file}"

python ${bad_mut_script} predict \
    -c ${config_file} \
    -f ${fasta_file} \
    -a ${msa_fasta} \
    -r ${msa_tree} \
    -s ${subs_file} \
    -o ${out_dir}/${subdir_name} \
    1> ${out_dir}/all_predict_log/${subdir_name}/${msa_name}_predict.log

if [[ -f ${out_dir}/${subdir_name}/${msa_name}_Predictions.txt ]]
then
    # Predict output file got written
    echo "Predict output file was written: ${out_dir}/${subdir_name}/${msa_name}_Predictions.txt"
    echo "Checking predict output file..."
    # More in depth check of predict file in cases where file was written but error still occured.
    #   e.g., "Error: HyPhy killed by signal 15", "Function call stack", etc.
    if grep -wq "Error" ${out_dir}/${subdir_name}/${msa_name}_Predictions.txt
    then
        # Move problematic prediction to subdirectory for easy troubleshooting
        mv ${out_dir}/${subdir_name}/${msa_name}_Predictions.txt ${out_dir}/${subdir_name}/problematic_predictions
        echo "The following prediction resulted in an error: ${out_dir}/${subdir_name}/problematic_predictions/${msa_name}_Predictions.txt"
        echo "Please troubleshoot."
    fi
else
    # Predict output file failed to get written to file, exiting
    echo "The following predict output file didn't successfully get written to a file: ${out_dir}/${subdir_name}/${msa_name}_Predictions.txt"
    echo "Exiting..."
    exit 2
fi





for i in $(seq 0 ${#MSA_FASTA_ARR[@]}); do if [[ "${MSA_FASTA_ARR[$i]}" == *"HORVU.MOREX.r3.2HG0128880.1"* ]]; then echo "$i : ${MSA_FASTA_ARR[$i]}"; fi; done

# HORVU.MOREX.r3.2HG0128600.1 # 0
HORVU.MOREX.r3.2HG0128750.1 # 1
HORVU.MOREX.r3.2HG0128800.1 # 2
HORVU.MOREX.r3.2HG0128880.1 # 3

# HORVU.MOREX.r3.2HG0124660.1 # 0
# HORVU.MOREX.r3.2HG0124740.1 # 1
# HORVU.MOREX.r3.2HG0124800.1 # 2
# HORVU.MOREX.r3.2HG0124810.1 # 3

#HORVU.MOREX.r3.1HG0074120.1 # 2
#HORVU.MOREX.r3.1HG0074020.1 # 0
#HORVU.MOREX.r3.1HG0074060.1 # 1
#HORVU.MOREX.r3.1HG0074140.1 # 3
#HORVU.MOREX.r3.1HG0074200.1 # 4
#HORVU.MOREX.r3.1HG0074210.1 # 5
