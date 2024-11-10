#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 20:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J Pairfq
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-53/26_S3_L001_R1_001/03_pairfq/BarLib_6.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-53/26_S3_L001_R1_001/03_pairfq/BarLib_6.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-53/26_S3_L001_R1_001/03_pairfq

# Base directory
base_dir="/home/hooimin/lu2024-17-19/RAAV-53"

# List of directories ending with _AAV
directories=($(find "$base_dir" -type d -name "*_R1_001"))

# Run pairfq on each directory
for dir in "${directories[@]}"; do
    for file1 in "$dir"/02_BarExt/*_6.fastq.gz; do
        # Input Files 1
        in_name_R1="$file1" 
        for file2 in "$dir"/04_library/*_LibExt_6.fastq.gz; do
        # Input Files 2
            in_name_R1L="$file2" 

            name_out="6_BarLib_paired_reads" # remember to change this name

            # Create temporary files
            out_name_R1=$(mktemp -p "$dir/03_pairfq/" "p006R1_pair_XXXXX.fastq.gz")
            out_name_R1L=$(mktemp -p "$dir/03_pairfq/" "p006R1L_pair_XXXXX.fastq.gz")
            out_name_R1_singlet=$(mktemp -p "$dir/03_pairfq/" "p006R1_singlet_XXXXX.fastq.gz")
            out_name_R1L_singlet=$(mktemp -p "$dir/03_pairfq/" "p006R1L_singlet_XXXXX.fastq.gz")

            # Run pairfq
            ~/perl5/bin/pairfq makepairs -c gzip -f $in_name_R1 -r $in_name_R1L -fp $out_name_R1 -rp $out_name_R1L -fs $out_name_R1_singlet -rs $out_name_R1L_singlet --stats 2>&1

            # Move and rename output files
            mv $out_name_R1 "$dir/03_pairfq/barcodes_$name_out.fastq.gz"
            mv $out_name_R1L "$dir/03_pairfq/library_$name_out.fastq.gz"

            # Cleanup temporary files
            rm $out_name_R1_singlet $out_name_R1L_singlet

        done
    done 
done   

