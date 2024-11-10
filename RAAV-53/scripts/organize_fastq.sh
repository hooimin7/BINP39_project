#!/bin/bash

# Base directory
base_dir="/home/hooimin/lu2024-17-19/RAAV-53"

# Find all .fastq.gz files in the base directory
fastq_files=($(find "$base_dir" -name "*fastq.gz"))

# Loop through each .fastq.gz file 
  for fastq_file in "${fastq_files[@]}"; do
    # Check if there are any .fastq.gz files
    [ -e "$fastq_file" ] || continue
    
    # Get the base name of the file (without the extension)
    base_name=$(basename "$fastq_file" .fastq.gz)
    
    # Create a directory with the same name as the base name
    mkdir -p "$base_dir/$base_name"
    
    # Move the .fastq file into the corresponding directory
    mv "$fastq_file" "$base_dir/$base_name/"
  done

echo "Folders created and files moved successfully."