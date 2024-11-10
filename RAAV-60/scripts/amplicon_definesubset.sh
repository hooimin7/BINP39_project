#!/bin/bash

#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 10:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J amplicon_subset
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-60/p006/fastq.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-60/p006/fastq.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-60/p006


# write this script to stdout-file 
cat $0

# Add seqtk to PATH
~/seqtk/seqtk


# Selection of real amplicons
echo "Selecting real amplicons......."

## Input Files
in_name_R1="p006_S3_L001_R1_001.fastq.gz" 
in_name_R2="p006_S3_L001_R2_001.fastq.gz" 

# This section searches the sequencing file and only select the files with valid amplicons
out_name_R1=$(mktemp -p . "p006R1_XXXXXXXXXX.fastq")
out_name_R2=$(mktemp -p . "p006R2_XXXXXXXXXX.fastq")

command_args="overwrite=true k=15 rcomp=f skipr2=t qhdist=0 maskmiddle=f hammingdistance=2 findbestmatch=f ordered=t threads=$(nproc) in=$in_name_R1 in2=$in_name_R2 outm=$out_name_R1 outm2=$out_name_R2 fliteral=GCCATCCTCTTATCTCGTGG"

~/bbmap/bbduk2.sh $command_args

# After bbduk2.sh command, compress the output files
gzip $out_name_R1
gzip $out_name_R2

# Update in_name variables to point to the compressed files
in_name_R1=$out_name_R1.gz
in_name_R2=$out_name_R2.gz

# Extraction of a subset
echo "Extracting subset......."

if [ "$run_subset" = true ]; then

    # Define the number of sequences to include in the subset
    subset_count=10000
    start_sequence=10001
    end_sequence=$((start_sequence + subset_count - 1))

    # List of input files
    input_files=("$in_name_R1" "$in_name_R2")

    for in_name in "${input_files[@]}"; do
        # seqtk or similar to sample a subset of fastq files
        # Define the output file name for the subset file
        out_name_sub="${in_name%.fastq.gz}_sub.fastq.gz"

        # Convert the input FASTQ file to a list of sequence names
        seqtk seq -A $in_name | awk '{if(NR%4==1) print substr($0,2)}' > seq_names.txt

        # Extract the sequence names for the subset
        sed -n "${start_sequence},${end_sequence}p" seq_names.txt > subset_seq_names.txt

        # Use seqtk subseq to create a subset of the input file
        seqtk subseq $in_name subset_seq_names.txt > $out_name_sub

        # Clean up the list of sequence names
        rm seq_names.txt subset_seq_names.txt
    done
fi

output_Reads1=$(($(gunzip -c $in_name_R1 | wc -l)/4))
output_Reads2=$(($(gunzip -c $in_name_R2 | wc -l)/4))

echo "Utilized sequences for R1: $output_Reads1"
echo "Utilized sequences for R2: $output_Reads2"