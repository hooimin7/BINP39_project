#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 150:00:00
#SBATCH --mem=190000
#SBATCH -N 2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J p006RNA_count
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-30/p006_AAV_02/06_traceBack/02_analysis/RNA.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-30/p006_AAV_02/06_traceBack/02_analysis/RNA.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-30/p006_AAV_02/06_traceBack

# write this script to stdout-file 
cat $0

# Load modules
module load GCC/12.3.0
module load R/4.3.2
R -e ".libPaths(c('~/MyRextensions', .libPaths()))"
export R_LIBS="~/MyRextensions:${R_LIBS}"

# ln -s ../../../RAAV-60/p005/06_library/fragments_5_FragBar_paired_again.fastq.gz fragments_5_FragBar_paired_again.fastq.gz
# ln -s ../../../RAAV-60/p005/06_library/barcodes_5_FragBar_paired_again.fastq.gz barcodes_5_FragBar_paired_again.fastq.gz


# Run R script to get the indexed FragmentsFinal in csv
Rscript ../../scripts/RNA_count_p006_INDchunk.R

# sort -t, -k1,1 02_analysis/found.p006_AAV_02.csv > 02_analysis/sorted_found.p006_AAV_02.csv
# awk -F, '!seen[$1]++' 02_analysis/sorted_found.p006_AAV_02.csv > 02_analysis/unique_sorted_found.p006_AAV_02.csv





