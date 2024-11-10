#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 10:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J collectRDS
#SBATCH -o /home/hooimin/lu2024-17-19/Data_RAAV/p005_data/03_normalize/allsamples.out
#SBATCH -e /home/hooimin/lu2024-17-19/Data_RAAV/p005_data/03_normalize/allsamples.err
#SBATCH -D /home/hooimin/lu2024-17-19/Data_RAAV/p005_data

# write this script to stdout-file 
cat $0

# Load modules
module load GCC/12.3.0
module load R/4.3.2
R -e ".libPaths(c('~/MyRextensions', .libPaths()))"
export R_LIBS="~/MyRextensions:${R_LIBS}"



# Run R script to get the indexed FragmentsFinal in csv
Rscript ../scripts/collectRDS_normalize_concatenatep005.R






