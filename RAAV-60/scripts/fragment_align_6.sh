#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 10:00:00
#SBATCH --mem=165G
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J fragment_align
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-60/p006/04_blast/02_analysis/frag_align.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-60/p006/04_blast/02_analysis/frag_align.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-60/p006/04_blast

# write this script to stdout-file 
cat $0

# Load modules
module load GCC/12.3.0
module load R/4.3.2
R -e ".libPaths(c('~/MyRextensions', .libPaths()))"
export R_LIBS="~/MyRextensions:${R_LIBS}"


# Run R script to get the indexed FragmentsFinal in csv
Rscript ../../scripts/Fragment_align_6.R





