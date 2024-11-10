#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH --mem=32G
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J bowtie2
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-60/p005/05_bowtie/bowtie.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-60/p005/05_bowtie/bowtie.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-60/p005/05_bowtie

# write this script to stdout-file 
cat $0

# Load modules
module load GCC/11.3.0
module load Bowtie2

# # Define the source directory where the original files are located
# source_dir="/home/hooimin/lu2024-17-19/RAAV-60/p005/04_blast"

# # Define the target directory (current directory in this case)
# target_dir=$(/home/hooimin/lu2024-17-19/RAAV-60/p005/05_bowtie)

# # Assuming you know the range of i or it's defined elsewhere in your script
# for i in {1..9}
# do
#     original_file="${source_dir}/subset${i}_5a.fna"
#     link_name="${target_dir}/subset${i}_5a.fna"
    
#     # Create the symbolic link
#     ln -s "$original_file" "$link_name"
# done

bowtie2-build ../04_blast/LUTdna.fna LUTdna_bowtie_index

  #   # Run bowtie2 for each subset (this works too, results in quite similar output)
  #   bowtie2 -N 1 -L 21 --gbar 21 --mp 4,2 -p 16 -x LUTdna_bowtie_index -U ../04_blast/fragments_5a_paired_unique.fastq.gz \
  # --un unaligned.fastq \
  # --al aligned.fastq \
  # --un-conc unaligned_conc.fastq \
  # --al-conc aligned_conc.fastq

bowtie2 -N 1 -L 21 --gbar 21 --mp 4,2 --very-sensitive -a -p 16 -x LUTdna_bowtie_index -U ../04_blast/fragments_5a_paired_unique.fastq.gz \
  --un unaligned.fastq \
  --al aligned.fastq \
  --un-conc unaligned_conc.fastq \
  --al-conc aligned_conc.fastq