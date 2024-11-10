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
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-60/p006/03_blast/blast.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-60/p006/03_blast/blast.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-60/p006/03_blast

# write this script to stdout-file 
cat $0

# Load modules
module load GCCcore/13.2.0
module load Python/3.11.5
module load Anaconda3
source config_conda.sh

conda activate seqkit

seqkit rmdup -s -i -o fragments_paired_unique.fastq.gz ../02_pairfq/fragments_paired_reads.fastq.gz

# seqkit seq -a fragments_paired_unique.fastq.gz > fragments_paired_unique.fasta

seqkit fq2fa fragments_paired_unique.fastq.gz > fragments_paired_unique.fasta

conda deactivate

# Run Python script
python ../../scripts/csv_to_fasta.py

conda activate blast

# Make database
makeblastdb -in LUTdna.fa -out LUTblast.db -dbtype nucl -parse_seqids -title LUT

# Run blastn with custom output format
# blastn -query fragments_paired_unique.fasta -db LUTblast.db -out fragments_paired_unique_blastn.out -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
# blastn -query fragments_paired_unique.fasta -db LUTblast.db -out fragments_paired_unique_blastn.out -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -evalue 1e-3 -word_size 7 -gapopen 2 -gapextend 1 -reward 2 -penalty -1
blastn -query fragments_paired_unique.fasta -db LUTblast.db -out fragments_paired_unique_blastn.out

