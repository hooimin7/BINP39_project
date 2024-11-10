#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 12:00:00
#SBATCH -N 2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J fragment_StarBar
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-72/S26_S1_L001_R1_001/06_traceBack/02_analysis/RNA_6.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-72/S26_S1_L001_R1_001/06_traceBack/02_analysis/RNA_6.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-72/

# write this script to stdout-file 
cat $0

# Load modules
module load GCC/12.3.0
module load R/4.3.2
R -e ".libPaths(c('~/MyRextensions', .libPaths()))"
export R_LIBS="~/MyRextensions:${R_LIBS}"

# Run R script 
Rscript scripts/RNA_count_p006_sub.R





