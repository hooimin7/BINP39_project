#!/bin/bash

# Base directory
base_dir="/home/hooimin/lu2024-17-19/RAAV-72"

# List of directories ending with _AAV
directories=($(find "$base_dir" -type d -name "*_R1_001"))

# List of folders to create
folders=("06_traceBack" "06_traceBack/03_p005_found" "06_traceBack/04_p007_found")

# Loop through each directory and create the folders
for dir in "${directories[@]}"; do
  for folder in "${folders[@]}"; do
    mkdir -p "$dir/$folder"
  done
done

echo "Folders created successfully."
