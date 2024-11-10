#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J fastqc
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-72/S26_S1_L001_R1_001/01_qualityTrim/01_fastqc/fastqc.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-72/S26_S1_L001_R1_001/01_qualityTrim/01_fastqc/fastqc.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-72/

# Base directory
base_dir="/home/hooimin/lu2024-17-19/RAAV-72"
module load FastQC/0.11.9-Java-11

# List of directories ending with _AAV
directories=($(find "$base_dir" -type d -name "*_R1_001"))

# Run FastQC on each directory
for dir in "${directories[@]}"; do
    echo "Running FastQC on $dir ..."
    fastqc --threads 6 --outdir "$dir/01_qualityTrim/01_fastqc" "$dir"/*.fastq.gz
done

# Write this script to stdout-file 
cat $0
