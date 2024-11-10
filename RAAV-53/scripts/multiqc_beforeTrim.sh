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
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-53/00_multiqc_result/01_beforeTrim/multiqcBT.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-53/00_multiqc_result/01_beforeTrim/multiqcBT.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-53/


# Base directory
base_dir="/home/hooimin/lu2024-17-19/RAAV-53"

# Load modules
module load GCC/12.2.0 OpenMPI/4.1.4
module load MultiQC/1.14

# Create a directory to store all fastqc.html files from 01_fastqc
multiqc_dir="$base_dir/00_multiqc_result/01_beforeTrim"
mkdir -p "$multiqc_dir"

# Collect all fastqc.zip files from 01_fastqc directories
find "$base_dir" -type f -path "*/01_fastqc/*_fastqc.zip" -exec cp {} "$multiqc_dir" \;

# Run MultiQC on the collected files
cd "$multiqc_dir"
multiqc .

# Write this script to stdout-file 
cat $0 