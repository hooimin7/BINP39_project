#!/bin/bash

#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J Starcode
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-60/p006/07_starcode/star.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-60/p006/07_starcode/star.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-60/p006/07_starcode

# Write this script to stdout-file 
cat $0

# Load modules
module load Anaconda3
source config_conda.sh
conda activate myenv

## Input Files
in_name_BC_paired="../06_library/barcodes_6_FragBar_paired_again.fastq.gz" # change input file name

# Barcode Outfiles
out_name_BC_star=$(mktemp -p . "BCsc_XXXXXX.txt")
num_cores=$(($(nproc) - 1))

# Decompress the input file
base_name=$(basename "$in_name_BC_paired" _FragBar_paired_again.fastq.gz)
processed_name1="${base_name}_reduced.fastq"
gunzip -c "$in_name_BC_paired" > "$processed_name1"

# Levenshtein distance of 1 (mismatch of 1 base pair), default is centered hamming distance
starcode -i "$processed_name1" -t $num_cores --print-clusters -d 1 -r5 -q -o "$out_name_BC_star"

# Ensure the data directory exists
mkdir -p 01_data

processed_name2="01_data/${base_name}_reduced.txt"
mv "$out_name_BC_star" "$processed_name2"

# Read the output file into a table
table_BC_sc=$(awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3}' "$processed_name2")

# Process the table
table_BC_sc=$(echo "$table_BC_sc" | awk -F"\t" '{split($2, a, ","); for (i in a) print $1, a[i]}')

# Calculate dropped BCs
reads_BC="../06_library/barcodes_6_FragBar_paired_again.fastq.gz"
# unique_reads_BC=$(zcat "$reads_BC" | grep -c "^@")

# Extract unique read barcodes, sort them, and count unique ones
unique_reads_BC=$(zcat "$reads_BC" | awk 'NR % 4 == 2' | sort | uniq | wc -l)
unique_table_BC_sc=$(echo "$table_BC_sc" | cut -f1 | sort | uniq | wc -l)
SC_droppedBC=$((unique_reads_BC - unique_table_BC_sc))

# Print the result
echo "Dropped BCs in Starcode: $SC_droppedBC"

# Save the result
echo "$table_BC_sc" > 01_data/scBC_DNA_pscAAVlib.txt

# Clean up
rm "$processed_name1"

