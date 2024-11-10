#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J fastqc
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-60/p007/01_qualityTrim/01_fastqc/fastqc.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-60/p007/01_qualityTrim/01_fastqc/fastqc.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-60/p007

# write this script to stdout-file 
cat $0

module load FastQC/0.11.9-Java-11

for i in /home/hooimin/lu2024-17-19/RAAV-60/p007/*.fastq.gz
do
    echo "Running $i ..."
    fastqc --threads 6 --outdir /home/hooimin/lu2024-17-19/RAAV-60/p007/01_qualityTrim/01_fastqc "$i"
done