#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 55:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J LibraryExt
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-53/26_S3_L001_R1_001/04_library/lib_7.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-53/26_S3_L001_R1_001/04_library/lib_7.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-53/26_S3_L001_R1_001/04_library

# write this script to stdout-file 
cat $0

# Load modules
module load Anaconda3
source config_conda.sh

# Base directory
base_dir="/home/hooimin/lu2024-17-19/RAAV-53"
num_cores=$(nproc --all)


# List of directories ending with _AAV
directories=($(find "$base_dir" -type d -name "*_R1_001"))

# Run bbduk2.sh on each directory
for dir in "${directories[@]}"; do
    for file in "$dir"/01_qualityTrim/*_subclean.fastq.gz; do
        # Input Files
        in_name_R1L="$file" 
   
        # Temporary output file
        out_name_R1L=$(mktemp -p . "p007R1L_lib_XXXXXXXXXX.fastq.gz")

        # Base name for processed file 
        base_name1=$(basename "$in_name_R1L" _subclean.fastq.gz)
        processed_name1="${base_name1}_LibExt_7.fastq.gz"

        # Run bbduk2.sh
        ~/bbmap/bbduk2.sh overwrite=true k=16 mink=13 hammingdistance=2 findbestmatch=t rcomp=f \
        qhdist=2 minavgquality=0 maxns=1 threads=$num_cores \
        in="$in_name_R1L" out="$out_name_R1L" lliteral="TGTGTCTATCGCAAGA" 
                                                        
        # Rename the output file
        mv "$out_name_R1L" "$dir/04_library/$processed_name1"

        # Run the Python script
        python /home/hooimin/lu2024-17-19/RAAV-53/scripts/libraryExt_p007.py "$dir/04_library/$processed_name1" "$dir/04_library/${base_name1}_53_007.fastq.gz"
    done
done





