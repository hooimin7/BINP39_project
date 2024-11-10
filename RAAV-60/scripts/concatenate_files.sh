#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 10:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J Concatenate
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-60/p007/02_fragBarExt/combine.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-60/p007/02_fragBarExt/combine.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-60/p007/02_fragBarExt

zcat p007a_S4_L001_R1_001_barExt_7.fastq.gz p007b_S5_L001_R1_001_barExt_7.fastq.gz | gzip > combined_p007_L001_R1_001_barExt.fastq.gz

zcat p007a_S4_L001_R2_001_fragExt_7.fastq.gz p007b_S5_L001_R2_001_fragExt_7.fastq.gz | gzip > combined_p007_L001_R2_001_fragExt.fastq.gz