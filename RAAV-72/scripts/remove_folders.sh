#!/bin/bash

# Base directory
base_dir="/home/hooimin/lu2024-17-19/RAAV-72"

# To remove a folder
# # List of directories ending with _AAV
# directories=($(find "$base_dir" -type d -name "*_R1_001"))

# # List of folders to be removed
# folders=("01_qualityTrim/01_fastqc/multiqc_01")

# # Loop through each directory and create the folders
# for dir in "${directories[@]}"; do
#   for folder in "${folders[@]}"; do
#     rm -rf "$dir/$folder"
#   done
# done

# To remove all files in a folder
# List of directories ending with _R1_001
directories=($(find "$base_dir" -type d -name "*_R1_001"))

# Loop through each directory and remove the contents of 02_fastqcAfterTrim
for dir in "${directories[@]}"; do
  rm -rf "$dir/06_traceBack/02_analysis/"*
done

echo "Folders removed successfully."
