#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH --mem=32G
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J blast
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-60/p005/04_blast/01_blastout/blastrep.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-60/p005/04_blast/01_blastout/blastrep.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-60/p005/04_blast

# write this script to stdout-file 
cat $0

# Load modules
module load Anaconda3
source config_conda.sh

conda activate blast


# Define the total number of lines and lines per file
total_lines=$(wc -l < fragments_5a_paired_unique.fna)
lines_per_file=$((total_lines / 9))

# Define the start and end line for the first file
start_line=1
end_line=$lines_per_file

for i in {1..9}
do
    # Use sed to extract lines from start_line to end_line
    sed -n "${start_line},${end_line}p" fragments_5a_paired_unique.fna > "subset${i}_5a.fna"

    # Run blastn for each subset
    # blastn -query "subset${i}_5a.fna" -db LUTblastNl -out "01_blastout/fragments_5a_paired_unique_blastn${i}.out" -word_size 21 -ungapped -num_threads 2
    blastn -query "subset${i}_5a.fna" -db LUTblastNl -out "01_blastout/fragments_5a_paired_unique_blastn${i}.out" -word_size 21 -ungapped -num_threads 2 -perc_identity 95.24
    
    # Calculate the total 'Query=' and 'No hits'
    total_query=$(grep -c 'Query=' "01_blastout/fragments_5a_paired_unique_blastn${i}.out")
    total_no_hits=$(grep -c 'No hits' "01_blastout/fragments_5a_paired_unique_blastn${i}.out")

    # Calculate the difference and print it
    difference=$(echo "$total_query - $total_no_hits" | bc -l)
    echo "Difference for subset${i}: $difference"

    # Update start_line and end_line for next iteration
    start_line=$((end_line + 1))
    end_line=$((start_line + lines_per_file - 1))

    # For the last file, make sure we go to the end of the file
    if [ $i -eq 9 ]; then
        end_line=$total_lines
    fi
done