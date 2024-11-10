#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH --mem=190000
#SBATCH -t 168:00:00
#SBATCH -N 2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J Motif
#SBATCH -o /home/hooimin/lu2024-17-19/Data_RAAV/p007_data/06_motif_72/hmm72.out
#SBATCH -e /home/hooimin/lu2024-17-19/Data_RAAV/p007_data/06_motif_72/hmm72.err
#SBATCH -D /home/hooimin/lu2024-17-19/Data_RAAV/p007_data

# # write this script to stdout-file 
# cat $0

# Load modules
module load GCCcore/13.2.0
module load Python/3.11.5
module load Anaconda3
source config_conda.sh
module load GCC/12.3.0
module load R/4.3.2
R -e ".libPaths(c('~/MyRextensions', .libPaths()))"
export R_LIBS="~/MyRextensions:${R_LIBS}"

ulimit -s 16384

# Run R script to get the indexed FragmentsFinal in csv
Rscript ../scripts/Peptide_cluster_p007_72.R






