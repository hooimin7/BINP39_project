#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --mem=128G
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J blast
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-60/p007/04_blast/01_blastout/blast_test.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-60/p007/04_blast/01_blastout/blast_test.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-60/p007/04_blast

# write this script to stdout-file 
cat $0

# Load modules
module load GCCcore/12.2.0
module load GCCcore/12.3.0
module load parallel/20230722
module load Anaconda3
source config_conda.sh

conda activate blast


# Define variables
fragments_unique_fa="fragments_7_paired_unique.fna"
blast_db="LUTblastNl"
blast_out_prefix="blast_output"

# Use parallel to process each subset
cat "$fragments_unique_fa" | parallel --block $(($(wc -l < $fragments_unique_fa)/$(nproc))) --recstart '>' --pipe \
    "blastn -max_target_seqs 25 -word_size 11 -num_threads 1 -outfmt 10 -db $blast_db -query - > ${blast_out_prefix}.{#}.csv"

# Concatenate all output files into a single file
cat ${blast_out_prefix}.*.csv > ${blast_out_prefix}.csv

# Clean up individual output files
rm ${blast_out_prefix}.*.csv

# Compress the BLAST output
gzip -c ${blast_out_prefix}.csv > 01_blastout/blast_output_7.csv.gz