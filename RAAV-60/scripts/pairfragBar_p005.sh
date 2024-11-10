#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 10:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J Pairfq
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-60/p005/03_pairfq/pairfq.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-60/p005/03_pairfq/pairfq.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-60/p005/03_pairfq

# Define input and output file names
infile_R1="../02_fragBarExt/combined_p005_L001_R1_001_barExt.fastq.gz"
infile_R2="../02_fragBarExt/combined_p005_L001_R2_001_fragExt.fastq.gz"
in_name_R1="R1_5_barEx.fastq"
in_name_R2="R2_5_fragEx.fastq"
name_out="5_paired_reads" # remember to change this name

# Gunzip infile_R1 and infile_R2 before proceeding
gunzip -c $infile_R1 > $in_name_R1
gunzip -c $infile_R2 > $in_name_R2

# Create temporary files
out_name_R1=$(mktemp -p . "p005R1_pair_XXXXX.fastq.gz")
out_name_R2=$(mktemp -p . "p005R2_pair_XXXXX.fastq.gz")
out_name_R1_singlet=$(mktemp -p . "p005R1_singlet_XXXXX.fastq.gz")
out_name_R2_singlet=$(mktemp -p . "p005R2_singlet_XXXXX.fastq.gz")

# Run pairfq
~/perl5/bin/pairfq makepairs -c gzip -f $in_name_R1 -r $in_name_R2 -fp $out_name_R1 -rp $out_name_R2 -fs $out_name_R1_singlet -rs $out_name_R2_singlet --stats 2>&1

# Move and rename output files
mv $out_name_R1 ./barcodes_$name_out.fastq.gz
mv $out_name_R2 ./fragments_$name_out.fastq.gz

# Cleanup temporary files
rm $out_name_R1_singlet $out_name_R2_singlet

