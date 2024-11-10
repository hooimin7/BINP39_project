#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 10:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J PairfqBr2
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-60/p005/06_library/fragBar2.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-60/p005/06_library/fragBar2.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-60/p005/06_library/

# Define input and output file names
in_name_R1P="barcodes_5_BarLib_paired_reads.fastq.gz"
in_name_R2P="../03_pairfq/fragments_5_paired_reads.fastq.gz"
name_out="5_FragBar_paired_again" # remember to change this name

# Create temporary files
out_name_R1P=$(mktemp -p . "p005R1P_pair_XXXXX.fastq.gz")
out_name_R2P=$(mktemp -p . "p005R2_pair_XXXXX.fastq.gz")
out_name_R1P_singlet=$(mktemp -p . "p005R1P_singlet_XXXXX.fastq.gz")
out_name_R2P_singlet=$(mktemp -p . "p005R2P_singlet_XXXXX.fastq.gz")

# Run pairfq
~/perl5/bin/pairfq makepairs -c gzip -f $in_name_R1P -r $in_name_R2P -fp $out_name_R1P -rp $out_name_R2P -fs $out_name_R1P_singlet -rs $out_name_R2P_singlet --stats 2>&1

# Move and rename output files
mv $out_name_R1P ./barcodes_$name_out.fastq.gz
mv $out_name_R2P ./fragments_$name_out.fastq.gz

# Cleanup temporary files
rm $out_name_R1P_singlet $out_name_R2P_singlet

