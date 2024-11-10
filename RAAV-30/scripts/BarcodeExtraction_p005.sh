#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J BarExt
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-30/p005_AAV_01/02_BarExt/extract_5.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-30/p005_AAV_01/02_BarExt/extract_5.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-30/p005_AAV_01/02_BarExt

# write this script to stdout-file 
cat $0

## Input Files
in_name_R1="../01_qualityTrim/AAV23-01_S1_R1_001_subclean.fastq.gz" # change input file name

# Barcode analysis
out_name_R1=$(mktemp -p . "p005R1_bar_XXXXXXXXXX.fastq.gz")
num_cores=$(nproc --all)

# Run bbduk2.sh
~/bbmap/bbduk2.sh overwrite=true k=16 mink=13 hammingdistance=2 findbestmatch=t rcomp=f \
qhdist=2 minavgquality=0 maxns=1 minlength=23 maxlength=24 threads=$num_cores \
in="$in_name_R1" out="$out_name_R1" lliteral="TGAACTTGGGACTTCG" rliteral="ATAACTTCGTATAATG" 2>&1

                                                                                                
base_name1=$(basename "$in_name_R1" _subclean.fastq.gz)
processed_name1="${base_name1}_barExt.fastq.gz"
mv "$out_name_R1" "$processed_name1"

