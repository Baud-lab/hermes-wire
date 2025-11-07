#!/bin/bash

non_host="/users/abaud/data/primary/P50/shallow"
dataset="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step1/dataset_raw"

# Change directory to the source folder
cd "$non_host"

# Loop through all files in the source folder
for file in *.fastq.gz; do
    # Check if the file name contains "BLANK", "_I1_", or "_I2_"
    if [[ "$file" == *BLANK* || "$file" == *_I1_* || "$file" == *_I2_* ]]; then
        continue  # Skip files containing "BLANK", "_I1_", or "_I2_"
    fi

    # Extract the first string before the separator "_"
    prefix=$(echo "$file" | cut -d "_" -f 1)

    # Check the length of the prefix and add "000" if it's less than 10 characters
    if [ ${#prefix} -lt 10 ]; then
        new_name="000$file"
    else
        new_name="$file"
    fi

    # Create a symbolic link in the destination folder
    ln -s "$non_host/$file" "$dataset/$new_name"
done

