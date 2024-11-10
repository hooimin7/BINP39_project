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
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-60/p006/02_pairfq/fastq.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-60/p006/02_pairfq/fastq.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-60/p006/02_pairfq

# Define input and output file names
in_name_R1="p006R1_uhkHLC8L1t_barExt.fastq"
in_name_R2="p006R2_HMKZQod827_fragExt.fastq"
name_out="paired_reads"

# Create temporary files
out_name_R1=$(mktemp "p006R1_pair_XXXXX.fastq.gz")
out_name_R2=$(mktemp "p006R2_pair_XXXXX.fastq.gz")
out_name_R1_singlet=$(mktemp "p006R1_singlet_XXXXX.fastq.gz")
out_name_R2_singlet=$(mktemp "p006R2_singlet_XXXXX.fastq.gz")

# Run pairfq
~/perl5/bin/pairfq makepairs -c gzip -f $in_name_R1 -r $in_name_R2 -fp $out_name_R1 -rp $out_name_R2 -fs $out_name_R1_singlet -rs $out_name_R2_singlet --stats 2>&1

# Move and rename output files
mv $out_name_R1 ./barcodes_$name_out.fastq.gz
mv $out_name_R2 ./fragments_$name_out.fastq.gz

# Cleanup temporary files
rm $out_name_R1_singlet $out_name_R2_singlet

echo "Total execution time:"