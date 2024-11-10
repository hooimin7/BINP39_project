#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J missing_check
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-60/p005/04_blast/02_analysis/missing.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-60/p005/04_blast/02_analysis/missing.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-60/p005/04_blast

# write this script to stdout-file 
cat $0

# Extract read IDs from the FASTQ file
zcat reorder_frag_5aFrgBr_pair_again.fastq_order.gz | awk 'NR % 4 == 1 {print substr($1, 2)}' > 02_analysis/fastq_read_ids.txt

# Extract read IDs from the BLAST output
zcat 01_blastout/blast_output.csv.gz | cut -d',' -f1 > 02_analysis/blast_read_ids.txt

# Sort the extracted read IDs
sort 02_analysis/fastq_read_ids.txt -o 02_analysis/fastq_read_ids.txt
sort 02_analysis/blast_read_ids.txt -o 02_analysis/blast_read_ids.txt

# Compare the read IDs and find missing ones
comm -23 02_analysis/blast_read_ids.txt 02_analysis/fastq_read_ids.txt > 02_analysis/missing_read_ids.txt