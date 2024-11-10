#!/bin/bash

# Base directory
base_dir="/home/hooimin/lu2024-17-19/RAAV-53"

# List of directories ending with _R1_001
directories=($(find "$base_dir" -type d -name "*_R1_001"))

# List of folders to delete files from
del_folders=("06_traceBack/02_analysis")

# Loop through each directory and create the folders
for dir in "${directories[@]}"; do
  for folder in "${del_folders[@]}"; do
    rm -r "$dir/$folder"/*
  done
done

echo "Files inside the specified folders deleted successfully."
