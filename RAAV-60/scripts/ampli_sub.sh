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
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-60/p005/fastq.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-60/p005/fastq.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-60/p005

# write this script to stdout-file 
cat $0

# Selection of real amplicons
# ============================
echo "Selecting real amplicons......."

## Input Files

in_name_R1="p006_S3_L001_R1_001.fastq.gz" 
in_name_R2="p006_S3_L001_R2_001.fastq.gz" 

# This section searches the sequencing file and only select the files with valid amplicons
out_name_R1=$(mktemp -p . "p006_XXXXXXXXXX.fastq.gz")
out_name_R2=$(mktemp -p . "p006_XXXXXXXXXX.fastq.gz")

command_args="overwrite=true k=15 rcomp=f skipr2=t qhdist=0 maskmiddle=f hammingdistance=2 findbestmatch=f ordered=t threads=$(nproc) in=$in_name_R1 in2=$in_name_R2 outm=$out_name_R1 outm2=$out_name_R2 fliteral=GCCATCCTCTTATCTCGTGG"

~/bbmap/bbduk2.sh $command_args

in_name_R1=$out_name_R1
in_name_R2=$out_name_R2

# Extraction of a subset
# ============================
echo "Extracting subset......."

function create_subset {
    local in_name_R1=$1
    local in_name_R2=$2
    local subset_count=100000

    local out_name_R1_sub="${in_name_R1%.fastq.gz}_sub.fastq.gz"
    local out_name_R2_sub="${in_name_R2%.fastq.gz}_sub.fastq.gz"

    seqtk sample -s100 $in_name_R1 $subset_count > $out_name_R1_sub
    seqtk sample -s100 $in_name_R2 $subset_count > $out_name_R2_sub

    # Check if the subset files were created
    if [[ -f $out_name_R1_sub && -f $out_name_R2_sub ]]; then
        echo "Subset files $out_name_R1_sub and $out_name_R2_sub created successfully."
    else
        echo "Error: Subset files not created."
    fi
}

if [ "$run_subset" = true ]; then
    create_subset "p006_S3_L001_R1_001.fastq.gz" "p006_S3_L001_R2_001.fastq.gz"
fi


output_Reads1=$(($(gunzip -c $in_name_R1 | wc -l)/4))
output_Reads2=$(($(gunzip -c $in_name_R2 | wc -l)/4))

echo "Utilized sequences for R1: $output_Reads1"
echo "Utilized sequences for R2: $output_Reads2"