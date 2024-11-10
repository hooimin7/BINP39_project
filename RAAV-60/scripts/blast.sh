#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J blast
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-60/p006/03_blast/fastq.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-60/p006/03_blast/fastq.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-60/p006/03_blast

# write this script to stdout-file 
cat $0

# load the modules required for running R
module load GCC/12.3.0
module load R/4.3.2
R -e ".libPaths(c('~/MyRextensions', .libPaths()))"
export R_LIBS="~/MyRextensions:${R_LIBS}"
module load Anaconda3
source config_conda.sh
conda activate blast
module load GCCcore/12.2.0
module load GCCcore/12.3.0
module load parallel/20230722
Rscript /home/hooimin/lu2024-17-19/RAAV-60/scripts/blast.R

