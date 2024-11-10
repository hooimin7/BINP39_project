#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J LibraryExt
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-30/p006_AAV_02/04_library/lib.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-30/p006_AAV_02/04_library/lib.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-30/p006_AAV_02/04_library

# write this script to stdout-file 
cat $0

# Input Files
in_name_R1L="../01_qualityTrim/AAV23-02_S2_R1_001_subclean.fastq.gz" 

# Library analysis
out_name_R1L=$(mktemp -p . "p006R1L_lib_XXXXXXXXXX.fastq.gz")
num_cores=$(nproc --all)

# Run bbduk2.sh

~/bbmap/bbduk2.sh overwrite=true k=16 mink=13 hammingdistance=2 findbestmatch=t rcomp=f \
qhdist=2 minavgquality=0 maxns=1 threads=$num_cores \
in="$in_name_R1L" out="$out_name_R1L" lliteral="TCCTCTTATCTCGTGG" 

base_name1=$(basename "$in_name_R1L" _subclean.fastq.gz)
processed_name1="${base_name1}_LibExt.fastq.gz"
mv "$out_name_R1L" "$processed_name1"

# Load modules
module load Anaconda3
source config_conda.sh

python /home/hooimin/lu2024-17-19/RAAV-30/scripts/libraryExt_p006.py

