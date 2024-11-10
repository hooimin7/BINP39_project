#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 168:00:00
#SBATCH --mem=190000
#SBATCH -N 2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J AllVsLib
#SBATCH -o /home/hooimin/lu2024-17-19/Data_RAAV/p007_data/04_plot/top25_AAV_lib_test/plotAAVLib_72.out
#SBATCH -e /home/hooimin/lu2024-17-19/Data_RAAV/p007_data/04_plot/top25_AAV_lib_test/plotAAVLib_72.err
#SBATCH -D /home/hooimin/lu2024-17-19/Data_RAAV/p007_data

# write this script to stdout-file 
cat $0

# Load modules
module load GCC/12.3.0
module load R/4.3.2
R -e ".libPaths(c('~/MyRextensions', .libPaths()))"
export R_LIBS="~/MyRextensions:${R_LIBS}"



# Run R script to get the indexed FragmentsFinal in csv
Rscript ../scripts/PlotAllGenesCoverage_copy2_ind_p007_AllVsLib.R






