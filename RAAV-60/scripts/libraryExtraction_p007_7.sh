#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH --mem=32G
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J LibraryExt_7
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-60/p007/06_library/lib_7b.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-60/p007/06_library/lib_7b.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-60/p007/06_library

# write this script to stdout-file 
cat $0

# Input Files
in_name_R1L="../01_qualityTrim/p007b_S5_L001_R1_001_subclean.fastq.gz" 

# Library analysis
out_name_R1L=$(mktemp -p . "p007R1L_lib_XXXXXXXXXX.fastq.gz")
num_cores=$(nproc --all)

# Run bbduk2.sh

~/bbmap/bbduk2.sh overwrite=true k=16 mink=13 hammingdistance=2 findbestmatch=t rcomp=f \
qhdist=2 minavgquality=0 maxns=1 threads=$num_cores \
in="$in_name_R1L" out="$out_name_R1L" lliteral="TGTGTCTATCGCAAGA" 
                                               

base_name1=$(basename "$in_name_R1L" _subclean.fastq.gz)
processed_name1="${base_name1}_LibExt_p007_7.fastq.gz"
mv "$out_name_R1L" "$processed_name1"

# Load modules
module load Anaconda3
source config_conda.sh

python /home/hooimin/lu2024-17-19/RAAV-60/scripts/libraryExt_p007_7.py

