#!/bin/bash
#SBATCH -A lu2024-2-56
#SBATCH -p lu32
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ho4588ta-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J QualityTrim_p006
#SBATCH -o /home/hooimin/lu2024-17-19/RAAV-60/p007/01_qualityTrim/qualityT.out
#SBATCH -e /home/hooimin/lu2024-17-19/RAAV-60/p007/01_qualityTrim/qualityT.err
#SBATCH -D /home/hooimin/lu2024-17-19/RAAV-60/p007/

# write this script to stdout-file 
cat $0

# Load modules
module load Anaconda3
source config_conda.sh
conda activate bbmap



for i in /home/hooimin/lu2024-17-19/RAAV-60/p007/*.fastq.gz
do
    echo "Running $i ..."
    base_name=$(basename "$i" .fastq.gz)
    output_file="/home/hooimin/lu2024-17-19/RAAV-60/p007/01_qualityTrim/${base_name}_subclean.fastq.gz"
    bbduk.sh -Xmx1g in="$i" out="$output_file" qtrim=rl trimq=15
done



