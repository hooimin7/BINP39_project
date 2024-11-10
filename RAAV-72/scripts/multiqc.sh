#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J multiqc
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-72/S26_S1_L001_R1_001/01_qualityTrim/01_fastqc/multiqc.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-72/S26_S1_L001_R1_001/01_qualityTrim/01_fastqc/multiqc.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-72/

# Base directory
base_dir="/home/hooimin/lu2024-17-19/RAAV-72"

# Load modules
module load GCC/12.2.0  OpenMPI/4.1.4
module load MultiQC/1.14

# List of directories ending with _R1_001
directories=($(find "$base_dir" -type d -name "*_R1_001"))

# Folder to create
folder="01_qualityTrim/01_fastqc/multiqc_01"

# Run MultiQC on each directory
for dir in "${directories[@]}"; do
    mkdir -p "$dir/$folder"
    cd "$dir/$folder"
    multiqc ../../01_fastqc/ .
done

# Write this script to stdout-file 
cat $0


