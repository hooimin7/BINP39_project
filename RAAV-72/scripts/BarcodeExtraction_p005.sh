#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 30:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J BarExt
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-72/S26_S1_L001_R1_001/02_BarExt/extract_5.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-72/S26_S1_L001_R1_001/02_BarExt/extract_5.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-72/S26_S1_L001_R1_001/02_BarExt

# write this script to stdout-file 
cat $0

# Base directory
base_dir="/home/hooimin/lu2024-17-19/RAAV-72"
num_cores=$(nproc --all)

# List of directories ending with _AAV
directories=($(find "$base_dir" -type d -name "*_R1_001"))

# Run barcode extraction on each directory
for dir in "${directories[@]}"; do
    for file in "$dir"/01_qualityTrim/*_subclean.fastq.gz; do
        # Input Files
        in_name_R1="$file" 
   
        # Temporary output file
        out_name_R1=$(mktemp -p . "p005R1_bar_XXXXXXXXXX.fastq.gz")

        # Base name for processed file 
        base_name1=$(basename "$in_name_R1" _subclean.fastq.gz)
        processed_name1="${base_name1}_barExt_5.fastq.gz"

        # Run bbduk2.sh
        ~/bbmap/bbduk2.sh overwrite=true k=16 mink=13 hammingdistance=2 findbestmatch=t rcomp=f \
        qhdist=2 minavgquality=0 maxns=1 minlength=23 maxlength=24 threads=$num_cores \
        in="$in_name_R1" out="$out_name_R1" lliteral="TGAACTTGGGACTTCG" rliteral="ATAACTTCGTATAATG" 2>&1


        # Rename the output file
        mv "$out_name_R1" "$dir/02_BarExt/$processed_name1"
    done
done    




                                                                        

