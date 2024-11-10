#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 10:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J FragExt
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-60/p006/fastq.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-60/p006/fastq.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-60/p006

# write this script to stdout-file 
cat $0

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

# Update in_name variables to point to the compressed files
in_name_R1=$out_name_R1
in_name_R2=$out_name_R2

# Fragment analysis
out_name_R2=$(mktemp -p . "p006R2_frg_XXXXXXXXXX.fastq.gz")
num_cores=$(nproc --all)

~/bbmap/bbduk2.sh overwrite=true k=12 mink=12 rcomp=f qhdist=1 maskmiddle=t \
hammingdistance=2 findbestmatch=t minlength=18 maxlength=24 ordered=t \
threads=$num_cores in="$in_name_R2" out="$out_name_R2" \
lliteral="CCAGAGAGGCAACGCT" rliteral="GCCAGACAAGCAGCTA" 2>&1

base_name2=$(basename "$in_name_R2" .fastq)
processed_name2="${base_name2}_fragExt.fastq"
mv "$out_name_R2" "$processed_name2"

# After bbduk2.sh command, compress the output files
gzip $processed_name2

# Barcode analysis
out_name_R1=$(mktemp -p . "p006R1_bar_XXXXXXXXXX.fastq.gz")
num_cores=$(nproc --all)

# Run bbduk2.sh
~/bbmap/bbduk2.sh overwrite=true k=12 mink=12 hammingdistance=2 findbestmatch=t rcomp=f \
qhdist=1 minavgquality=0 maxns=0 minlength=18 maxlength=24 threads=$num_cores \
in="$in_name_R1" out="$out_name_R1" lliteral="GTCTATCGCAAGACTA" rliteral="ATAACTTCGTATAATG" 2>&1

base_name1=$(basename "$in_name_R1" .fastq)
processed_name1="${base_name1}_barExt.fastq"
mv "$out_name_R1" "$processed_name1"

# After bbduk2.sh command, compress the output files
gzip $processed_name1