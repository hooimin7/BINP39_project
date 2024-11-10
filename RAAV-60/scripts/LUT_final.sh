#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 10:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J LUT
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-60/p006/06_library/01_LUT/LUT.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-60/p006/06_library/01_LUT/LUT.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-60/p006/06_library

# Define the names of the input files and the output file
barcode_file="barcodes_6_FragBar_paired_again.fastq.gz"
fragment_file="fragments_6_FragBar_paired_again.fastq.gz"
library_file="library_6_BarLib_paired_reads.fastq.gz"
output_file="01_LUT/6_fragment_barcode_library.txt"

# Extract the identifiers and sequences from the fragment file
zcat $fragment_file | awk 'NR%4==1 || NR%4==2' > fragment_ids_seqs.txt

# Extract the identifiers and sequences from the barcode file
zcat $barcode_file | awk 'NR%4==1 || NR%4==2' > barcode_ids_seqs.txt

# Extract the identifiers and sequences from the library file
zcat $library_file | awk 'NR%4==1 || NR%4==2' > library_ids_seqs.txt

# Interleave the identifiers and sequences within each file
awk '{getline x; print $0, x}' fragment_ids_seqs.txt > fragment_interleaved.txt
awk '{getline x; print $0, x}' barcode_ids_seqs.txt > barcode_interleaved.txt
awk '{getline x; print $0, x}' library_ids_seqs.txt > library_interleaved.txt

# Combine the three interleaved files into a single file with six columns
echo -e "Fragment_ID\tFragment_Seq\tBarcode_ID\tBarcode_Seq\tLibrary_ID\tLibrary_Seq" > $output_file
paste -d '\t' fragment_interleaved.txt barcode_interleaved.txt library_interleaved.txt >> $output_file

# Remove temporary files
rm fragment_ids_seqs.txt barcode_ids_seqs.txt library_ids_seqs.txt fragment_interleaved.txt barcode_interleaved.txt library_interleaved.txt