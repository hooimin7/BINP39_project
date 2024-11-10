#!/bin/bash

# write this script to stdout-file 
cat $0

# Load modules
module load GCC/12.3.0
module load R/4.3.2
R -e ".libPaths(c('~/MyRextensions', .libPaths()))"
export R_LIBS="~/MyRextensions:${R_LIBS}"
module load GCC/12.3.0
module load SAMtools/1.18
# Load conda #
module load Anaconda3
source config_conda.sh
conda activate bowtie2


# Run R script to get the indexed FragmentsFinal in csv
Rscript ../../scripts/S3_libraryIdentification_mod.R
