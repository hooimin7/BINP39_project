#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J QualityTrim
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-72/S26_S1_L001_R1_001/01_qualityTrim/qualityT.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-72/S26_S1_L001_R1_001/01_qualityTrim/qualityT.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-72/

# write this script to stdout-file 
cat $0

# Load modules
module load Anaconda3
source config_conda.sh
conda activate bbmap

# Base directory
base_dir="/home/hooimin/lu2024-17-19/RAAV-72"

# List of directories ending with _AAV
directories=($(find "$base_dir" -type d -name "*_R1_001"))

# Run quality trim 
for dir in "${directories[@]}"; do
    # Process each fastq.gz file in the directory
    for file in "$dir"/*.fastq.gz; do
     # Get the base name of the file (without the extension)
    base_name=$(basename "$file" .fastq.gz)
    output_file="$dir/01_qualityTrim/${base_name}_subclean.fastq.gz"
    # Run bbduk.sh
    bbduk.sh -Xmx1g in="$file" out="$output_file" qtrim=rl trimq=15 
    done
done


