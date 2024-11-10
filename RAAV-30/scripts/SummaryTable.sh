#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J SummaryT
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-30/p007_AAV_03/06_traceBack/01_table/summary.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-30/p007_AAV_03/06_traceBack/01_table/summary.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-30/p007_AAV_03/06_traceBack/

# write this script to stdout-file 
cat $0

# Load modules
module load Anaconda3
source config_conda.sh

python /home/hooimin/lu2024-17-19/RAAV-30/scripts/summaryTableplus.py

