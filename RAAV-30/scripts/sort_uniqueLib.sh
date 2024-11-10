#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J Sort_unique
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-30/p007_AAV_03/04_library/02_sortUniq/sortLib_7.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-30/p007_AAV_03/04_library/02_sortUniq/sortLib_7.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-30/p007_AAV_03/04_library

# Ensure the data directory exists
mkdir -p 02_sortUniq

# write this script to stdout-file 
cat $0

# Input FASTQ file
input_fastq="../03_pairfq/library_7_BarLib_paired_reads.fastq.gz"

# Output file for unique sorted sequences with counts
output_file="02_sortUniq/library_7_BarLib_paired_reads_sortU.txt"

# Extract, sort, and count unique sequences
zcat "$input_fastq" | awk 'NR % 4 == 2' | sort | uniq -c | sort -nr > "$output_file"

# Sum the counts
total_count=$(awk '{sum += $1} END {print sum}' "$output_file")

echo "Total count of unique sequences: $total_count"
