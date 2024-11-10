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
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-60/p007/04_blast/blast.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-60/p007/04_blast/blast.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-60/p007/04_blast

# write this script to stdout-file 
cat $0

# Load modules
module load GCCcore/13.2.0
module load Python/3.11.5
module load Anaconda3
source config_conda.sh
module load GCC/12.3.0
module load R/4.3.2
R -e ".libPaths(c('~/MyRextensions', .libPaths()))"
export R_LIBS="~/MyRextensions:${R_LIBS}"

conda activate seqkit

# Run extracted fragments to get the fasta format
seqkit rmdup -s -P -o fragments_7_paired_unique.fastq.gz ../06_library/fragments_7_FragBar_paired_again.fastq.gz

seqkit fq2fa fragments_7_paired_unique.fastq.gz > fragments_7_paired_unique.fna

conda deactivate

# Run R script to get the indexed FragmentsFinal in csv
Rscript ../../scripts/fragmentsFinal.R

# Run Python script
python ../../scripts/csv_to_fasta.py

conda activate blast

# Make database
makeblastdb -in LUTdna.fna -out LUTblastNl -dbtype nucl -parse_seqids -title LUT


