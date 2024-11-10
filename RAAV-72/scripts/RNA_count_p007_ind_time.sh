#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 10:00:00
#SBATCH --mem=190000
#SBATCH -N 2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J RNA_count
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-72/S59_S8_L001_R1_001/06_traceBack/04_p007_found/RNA_7.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-72/S59_S8_L001_R1_001/06_traceBack/04_p007_found/RNA_7.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-72/S59_S8_L001_R1_001

# write this script to stdout-file 
cat $0

# Load modules
module load GCC/12.3.0
module load R/4.3.2
R -e ".libPaths(c('~/MyRextensions', .libPaths()))"
export R_LIBS="~/MyRextensions:${R_LIBS}"

# Run R script to get the indexed FragmentsFinal in csv
Rscript ../scripts/RNA_count_p007_INDchunk_time.R 
