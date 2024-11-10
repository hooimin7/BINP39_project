#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 10:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J Concatenate_lib
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-60/p005/06_library/combine.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-60/p005/06_library/combine.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-60/p005/06_library

zcat library_005a.fastq.gz library_005b.fastq.gz | gzip > combined_p005_LibExt.fastq.gz

