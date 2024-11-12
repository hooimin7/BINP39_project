# BINP39_project
# Author's note #
Environment: 

#### Load conda 24.1.2 in Lunarc ####
```bash
module load Anaconda3
source config_conda.sh
```
#### trimmomatic version 0.39 #### 
```bash
conda install bioconda::trimmomatic
conda activate trimmomatic
```
#### FastQC: FastQC/0.11.9-Java-11 #### 
available on lunarc
```bash
module load FastQC/0.11.9-Java-11
```
#### MultiQC: MultiQC/1.14 #### 
(may affect bowtie2 and jellyfish) available on lunarc
```bash
module load MultiQC/1.14
```
#### bowtie2 version 2.5.4 ####
```bash
conda create -n myenv
conda install bioconda::bowtie2
```
#### bbmap version 38.4 ####
```bash
conda create -n bbmap
conda install bioconda::bbmap
```
#### Install local::lib and Tie::Hash::DBD ####
```bash
cpan local::lib
eval "$(perl -I$HOME/perl5/lib/perl5 -Mlocal::lib)"
cpan Tie::Hash::DBD
```
#### pairfq 0.17.0 ####
[pairfq](https://github.com/sestaton/Pairfq/releases) # downloaded version
put in the ./bin
```bash
unzip Pairfq-0.17.0.zip -d /home/hooimin/bin
# Install Pairfq to local::lib
cd ~/bin/Pairfq-0.17.0
perl Makefile.PL PREFIX=~/perl5 LIB=~/perl5/lib/perl5
perl Makefile.PL INSTALL_BASE=~/perl5
make
make test
make install PREFIX=~/perl5 LIB=~/perl5/lib/perl5
```
#### blastn 2.15.0+ ####
```bash
conda install -c bioconda blast 
```
#### seqkit version 2.8.2 ####
```bash
conda install bioconda::seqkit
```
#### HHMER 3.4 ####
```bash
conda install bioconda::hmmer
```
#### Hammock v_1.2.0 ####
```bash
wget https://github.com/krejciadam/hammock/releases/download/v1.2.0/Hammock_v_1.2.0.7z
# Extract the archive:
7z x ~/Hammock_v_1.2.0.7z -o~/Hammock_v_1.2.0
# List the contents of the extracted directory
ls -R ~/Hammock_v_1.2.0
```
#### Weblogo 0.0.0 ####
```bash
conda install bioconda::weblogo
```
#### Ghostscript 9.18 (2015-10-05) ####
```bash
conda install bioconda::ghostscript
```
#### starcode version 1.4 ####
```bash
conda create -n myenv
conda install bioconda::starcode
```
#### entrez-direct version 22.4 ####
```bash
conda install bioconda::entrez-direct
```
#### Install jq version 1.5 ####
```bash
conda install bioconda::jq
```
#### R environment ####
```bash
mkdir ~/MyRextensions # create a directory for R packages
> .libPaths(c('~/MyRextensions', .libPaths())) 
> install.packages('BiocManager') # install package inside my diskspace
> BiocManager::install("Biostrings") # example of installing packages
```
#### Load Python in Lunarc ####
```bash
module load GCCcore/13.2.0
module load Python/3.11.5
```
#### Load Biopython ####
```bash
module load GCC/12.2.0  OpenMPI/4.1.4
module load Biopython/1.81
```
#### parallel/20230722 ####
```bash
module load GCCcore/12.2.0
module load GCCcore/12.3.0
module load parallel/20230722
```
#### google cloud ####
```bash
curl https://sdk.cloud.google.com | bash
### path saved at .bashrc, reset the shell ###
exec -l $SHELL
```
#### Add PATH for the Google Cloud SDK ####
```bash
if [ -f '~/google-cloud-sdk/path.bash.inc' ]; then
    . '~/google-cloud-sdk/path.bash.inc'
fi
```
#### Add seqtk to PATH ####
```bash
export PATH=$PATH:~/seqtk
```
# 1. Raw data processing #
## RAAV-60: p005, apply for p006, p007 ##
Some with similar scripts, some with special scripts for each group 
Fastq files downloaded from Google Cloud were saved in RAAV-60: p005, p006 and p007
## fastqc before trim ##
Check the fastq files quality, p005 and p007 have duplication of files (two pairs pf forward, two pairs of reverse files), p006 has only one pair of forward and reverse files
```bash
sbatch ~/RAAV-60/scripts/fastqc.sh
```
## 01_qualiyTrim ##
```bash
bbduk.sh -Xmx1g in=p005a_S1_L001_R1_001.fastq.gz out=p005a_S1_L001_R1_001_subclean.fastq.gz qtrim=rl trimq=15 # for testing one #
sbatch ~/RAAV-60/scripts/qualityTrim.sh
```
## 01_qualityTrim/02_fastqcAfterTrim ## 
```bash
sbatch ~/RAAV-60/scripts/fastqcAfterTrim.sh
```
## MultiQc ##
```bash
module load GCC/12.2.0  OpenMPI/4.1.4
module load MultiQC/1.14
# stands at 01_fastqc/multiqc_01 #
multiqc ../../01_fastqc/ .
# stands at 01_fastqc/multiqc_01_aftertrim #
multiqc ../../02_fastqcAfterTrim/ .
```
# 2. Extraction of Barcodes, gene fragments and library IDs #
## RAAV-60: p005, 02_fragBarExt ##
#### scripts ####
```bash
sbatch ~/RAAV-60/scripts/fragmentExtraction_p005.sh
sbatch ~/RAAV-60/scripts/fragmentExtraction_p006_6.sh
sbatch ~/RAAV-60/scripts/fragmentExtraction_p007_7.sh
```
#### Barcode and fragment extraction #### 
It takes around 14 hours clock time to run, examples for p005
Fragment extarction:
```bash
~/bbmap/bbduk2.sh overwrite=true k=16 mink=13 hammingdistance=2 findbestmatch=t rcomp=f \
qhdist=2 minavgquality=0 maxns=1 maskmiddle=t minlength=21 maxlength=21 ordered=t threads=$num_cores \
in="$in_name_R2" out="$out_name_R2" lliteral="CCAGAGAGGCAACGCT" rliteral="GCCAGACAAGCAGCTA" 2>&1
```
Barcode extraction: -- use p005, p006, p007 own template for extraction
```bash
~/bbmap/bbduk2.sh overwrite=true k=16 mink=13 hammingdistance=2 findbestmatch=t rcomp=f \
qhdist=2 minavgquality=0 maxns=1 minlength=23 maxlength=24 threads=$num_cores \
in="$in_name_R1" out="$out_name_R1" lliteral="TGAACTTGGGACTTCG" rliteral="ATAACTTCGTATAATG" 2>&1
```

#### Concatenate extracted fragments and barcodes of p005a and p005b fastq files, of p007a and p007b fastq files ####
```bash
sbatch ~/RAAV-60/scripts/concatenate_files.sh 
zcat p005a_S1_L001_R1_001_barExt.fastq.gz p005b_S2_L001_R1_001_barExt.fastq.gz | gzip > combined_p005_L001_R1_001_barExt.fastq.gz
zcat p005a_S1_L001_R2_001_fragExt.fastq.gz p005b_S2_L001_R2_001_fragExt.fastq.gz | gzip > combined_p005_L001_R2_001_fragExt.fastq.gz
```

# 3. Pairfq pair matching #
## 03_pairfq ## 
Extracted fragment paired with extracted barcode, use pairfragBar_p006.sh for p006; use pairfragBar_p007.sh for p007 
The gunzip steps should be done for only once, one gzip step is removed in the script fragmentExtraction_p005.sh <br> 
***gunzip steps before pairfq*** <br>
The output files were zipped twice from bbduk2.sh 
```bash
sbatch ~/RAAV-60/scripts/pairfragBar_p005.sh 
```

Lookup table generation (For a general control of correct pairing of barcode to fragment), script example
```bash
sbatch ~/RAAV-60/scripts/LUT_p005.sh
```

# 4. Align against the NCBI genes using blast #
## 04_blast ##
Run mkblast with 06_library/fragments_5aFragBar_paired_again.fastq.gz ; Apply to p006 and p007
```bash
sbatch ~/RAAV-60/scripts/mkblast_HM.sh
```
There are two more scripts in the mkblast script, it needs fragmentsFinal.csv
```bash
Rscript ../../scripts/fragmentsFinal.R
python ../../scripts/csv_to_fasta.py
```
Running the blast for each subsets -- testing purposes
```bash
sbatch ~/RAAV-60/scripts/blast_HM_p005.sh 
```
blast_HM_p005.sh was modified to allow 20/21 similarity, but the same result was obtained

Blast in outfmt csv - for each group: p005, p006, p007
```bash
mkdir 01_blastout in 04_blast
sbatch ~/RAAV-60/scripts/blast_HM_outfmt_prefix.sh
```

# 5. Reverse mapping of CustumArray oligos to original proteins #
alignedLibraries.rda was created 
## 05_bowtie ##
```bash
sbatch ~/RAAV-60/scripts/S3_libraryIdentification.sh
Rscript ../../scripts/S3_libraryIdentification_mod.R
```
# 6. Library id pairing #
## 06_library ##
Use bbduk2 method
```bash
sbatch ~/RAAV-60/scripts/libraryExtraction_p005.sh
python ../../scripts/libraryExt_p006.py
```
Concatenate p005a and p005b library extracted, apply to p007a p007b
```bash
sbatch ~/RAAV-60/scripts/concatenate_libExtract.sh
```
Sort unique library (after paired), apply to all p005,p006, p007
***please mkdir 02_sortUniq before, otherwise it will fail***
Note: use input file ended with .fastq.gz (three nucleotides library id)
```bash
sbatch ~/RAAV-60/scripts/sort_uniqueLib.sh
```
Pairfq barcode and library p005a (Paired barcode paired with extracted library)
```bash
sbatch ~/RAAV-60/scripts/pairBarLib_p005.sh
in_name_R1="../03_pairfq/barcodes_5a_paired_reads.fastq.gz"
in_name_R1L="library_005a.fastq.gz" # (after combination: use combined_p005_LibExt.fastq.gz)
```
Note: for library p006: use library_006.fastq.gz as input file

Pairfq library-barcode with fragment p005a ### (Paired barcode-library paired with the paired fragment)
```bash
sbatch ~/RAAV-60/scripts/pairfragBar2_p005.sh
in_name_R1P="barcodes_5aBarLib_paired_reads.fastq.gz"
in_name_R2P="../03_pairfq/fragments_5a_paired_reads.fastq.gz"
```
Make a table of fragment, barcodes and library id
LUT frag-bar-lib, apply for p005, p006, p007
```bash
mkdir 01_LUT
sbatch ~/RAAV-60/scripts/LUT_final.sh
```
# 7. Starcode based barcode reduction #
## 07_starcode ## 
Apply for p005, p006, p007
```bash
sbatch ~/RAAV-60/scripts/starcode.sh
```

----------------------------------------------------------------
# 8. Raw data processing for RAAV-60: p006 #
fastqc before trim same as # 1
quality trim, fastqc after trim same as # 1
Start with library extraction, to check which library is involved, same as # 6
```bash
sbatch ~/RAAV-60/scripts/libraryExtraction_p006_6.sh 
```
Sort the unique library id from extracted library 
```bash
sbash ~/RAAV-60/scripts/sort_uniqueLib.sh
```
According to the sorting of unique: p006 library extracted showed 'TAC' is the highest, 
will use the sequence 'TCTTATCTCGTGGTAC' from RAAV-template sequence

Fragment extraction: with p006 barcode template sequence is correct
434710288 Aug 13 10:32 p006_S3_L001_R1_001_barExt_6.fastq.gz # (this is correct)
99287618 Aug 13 10:41 p006_S3_L001_R1_001_barExt_7.fastq.gz # (this is incorrect)

Since p006_S3_L001_R1_001_barExt_6.fastq.gz is correct, Proceed pairfq with p006_6 output

# 9. Raw data processing for RAAV-60: p007 #
p007 frag extraction,  use p007 barcode template sequence
```bash
sbatch ~/RAAV-60/scripts/fragmentExtraction_p007_7.sh
# Concatenate 7a and 7b (R1), 7a and 7b (R2)
sbatch ~/RAAV-60/scripts/concatenate_files.sh 
# p007 lib extraction #
sbatch ~/RAAV-60/scripts/libraryExtraction_p007_7.sh 
# Pairfq barcode and library p007a # (Paired barcode paired with extracted library)
sbatch ~/RAAV-60/scripts/pairBarLib_p007.sh
# Pairfq library-barcode with fragment p007a ### (Paired barcode-library paired with the paired fragment)
sbatch ~/RAAV-60/scripts/pairfragBar2_p007.sh
```
-------------------------------------------------------
# 10. Raw data processing #
This section will cover for RAAV-30, RAAV-53, RAAV-72 
## RAAV-30 virus ###
Create folders 
```bash
chmod +x create_folders.sh
./create_folders.sh 
```
## 01_qualityTrim ##
fastqc before trimming
```bash
sbatch ~/RAAV-30/scripts/fastqc_folders.sh
```
fastqc after trimming 
```bash
./create_folders.sh
02_fastqcAfterTrim
mkdir 02_fastqcAfterTrim 
sbatch ~/RAAV-30/scripts/fastqc_folders_afterTrim.sh
```
MultiQC (after trimmming), this applies to RAAV-30, RAAV-53, RAAV-72
```bash
sbatch ~/RAAV-30/scripts/multiqc_afterTrim.sh
```
Removed multiqc_01 that created earlier, this applies to RAAV-53, RAAV-72
```bash
./remove_folders.sh
```
## RAAV-53 and RAAV-72 tissue ##
Organize folder based on the fastq files name 
```bash
chmod +x create_folders.sh
./organize_fastq.sh
```
Create folders in each folder 
```bash
chmod +x create_folders.sh
./create_folders.sh
```
## 01_fastqc before trimming ##
```bash
sbatch ~/RAAV-53/scripts/fastqc_folders.sh
# MultiQC (before trimmming), this applies to RAAV-30, RAAV-53, RAAV-72
mkdir 00_multiqc_result
mkdir 01_beforeTrim 02_afterTrim
sbatch ~/RAAV-30/scripts/multiqc_beforeTrim.sh
# Quality trim, this applies to RAAV-30,  RAAV-53, RAAV-72
sbatch ~/RAAV-30/scripts/qualityTrim.sh 
# Trimmomatic run test comparison using RAAV-53
sbatch ~/RAAV-53/scripts/Trimmomatic.sh
```
***Due to a loss of unique reads in RAAV-53, RAAV-72 after Trimmomatic***
quality trim (bbduk) was used and applied to RAAV-53 and RAAV-72 

## 02_BarExt RAAV-30 ## 
```bash
sbatch ~/RAAV-30/scripts/BarcodeExtraction_p005.sh
sbatch ~/RAAV-30/scripts/BarcodeExtraction_p006_7.sh # (run again using this code due to too little output, use p007 template)
sbatch ~/RAAV-30/scripts/BarcodeExtraction_p006.sh 
sbatch ~/RAAV-30/scripts/BarcodeExtraction_p007_6.sh # (run again using this code due to too little output, use p006 template)
sbatch ~/RAAV-30/scripts/BarcodeExtraction_p007.sh
```
This library id extraction was done before 03_Pairfq
## 04_library RAAV-30 ## 
```bash
sbatch ~/RAAV-30/scripts/LibraryExtraction_p005.sh
sbatch ~/RAAV-30/scripts/LibraryExtraction_p006.sh
sbatch ~/RAAV-30/scripts/LibraryExtraction_p007.sh
```
To sort the unique library id
## 04_library/02_sortUniq/ RAAV-30 sort-uniq ##
```bash
sbatch ~/RAAV-30/scripts/sort_uniqueLib.sh 
```
## 03_Pairfq RAAV-30 barcode and library ## 
p005_AAV_01 (extracted barcode paired with extracted library)
```bash
sbatch ~/RAAV-30/scripts/pairBarLib_p005.sh
# p006_AAV_02 use the extarcted library from template p007, see ## 02_BarExt RAAV-30 ##
sbatch ~/RAAV-30/scripts/pairBarLib_p006_7.sh 
# p007_AAV_03 #### use the extarcted library from template p006, see ## 02_BarExt RAAV-30 ##
sbatch ~/RAAV-30/scripts/pairBarLib_p007_6.sh 
```

## 02_BarExt RAAV-53 ## 
Apply to RAAV-53, RAAV-72 
```bash
./create_folders.sh # Create: folders=("02_BarExt" "03_pairfq" "04_library" "05_starcode")
sbatch ~/RAAV-53/scripts/BarcodeExtraction_p005.sh 
```

## 04_library RAAV-53 ## 
(RAAV-53 has 20 samples, 3 long reads, the rest are short reads, this script applies to RAAV-72)
```bash
sbatch ~/RAAV-53/scripts/libraryExtraction_p005.sh 
```

## 03_Pairfq barcode and library RAAV-53, RAAV-72 ## 
(extracted barcode paired with extracted library)
Both RAAV-53 and RAAV-72 use the same script
```bash
sbatch ~/RAAV-53/scripts/pairBarLib_p005.sh
```

## 05_starcode ## 
Apply to RAAV-30, p005_AAV_01, p006_AAV_02, p007_AAV_03
```bash
sbatch ~/RAAV-30/scripts/starcode.sh
```
## 05_starcode ## 
Apply to RAAV-53, RAAV-72
```bash
sbatch ~/RAAV-53/scripts/starcode_5.sh 
```
----------------------------------------------------

# 11. Data analysis - found fragments #
## 11.1 p006 without LUTdna ##
The whole pipeline from creating the multipleContfragmentsCompleteNoLUT.rda to circular plot creation, skip matching with the fragmentsFinal.csv
```bash
sbatch ~/RAAV-60/scripts/combineStarBar_p006.sh
Rscript ../../scripts/combinefragStarBar_p006.R

sbatch ~/RAAV-53/scripts/RNA_count_p006noLUT.sh
Rscript scripts/RNA_count_p006_noLUT.R

sbatch ~/RAAV-72/scripts/RNA_count_p006noLUT.sh
Rscript scripts/RNA_count_p006_noLUT.R

Rscript scripts/cleanName.R # run at the terminal

sbatch ~/Data_RAAV/scripts/normalize_noLUT.sh
Rscript ../scripts/normalize_noLUT.R

sbatch ~/Data_RAAV/scripts/VennPlot.sh
Rscript ../scripts/VennPlot_p006_noLUT.R
```

## 11.2 Calculate consensus alignment of chimeric barcodes ##
Consists of multipleContfragmentsComplete.csv (Fragment+Barcode+LUT created with mCount, output.Table created for next step RNA tissue)
This applies to (RAAV-60: p005, p006, p007)
```bash
sbatch ~/RAAV-60/scripts/combineStarBarLUT_p005.sh
Rscript combinefragStarBarLUT_p005.R
# To calculate barcode, fragment frequency
sbatch ~/RAAV-60/scripts/fragBarFreq.sh  # (use this script to remove repetitve fragments)
```
## 11.3 Calculate alignment percentage ## 
RAAV-60: -p005, p006, p007 (take at least one hour to run for p005), p006 not aligned as it is not from fragmentsFinal.csv, scBC_DNA_pscAAVlib_p005.rda created
```bash
sbatch ~/RAAV-60/scripts/fragment_align.sh
Rscript Fragment_align.R

# Table creation using knitr
sbatch ~/RAAV-60/scripts/table_creation.sh 
```
***p007 was too big, have to create a smaller subset for the markdown***

## 11.4 RAAV-30 barcode analysis ##
```bash
sbatch ~/RAAV-30/scripts/combineStarBarLUT_p005.sh 

# Summary Table (RAAV-30, RAAV-53, RAAV-72)
# Table for capturing the statistic of different mrna (only done for RAAV-30) #
sbatch ~/RAAV-30/scripts/summaryTable.sh
python ~/RAAV-30/scripts/summaryTable.py
```

## 11.5 Generate a complete library range object ## 
To make a .rds file (combining multipleContfragmentsComplet with LUTdna) -- apply to (RAAV-60: p005, p007), a completeLibraryRanges.rds was created 
```bash
sbatch ~/RAAV-60/scripts/coverage_p005.sh 
sbatch ~/RAAV-60/scripts/coverage_p006.sh 
sbatch ~/RAAV-60/scripts/coverage_p007.sh 
```
*RNA_count (RAAV-30): p006_AAV_02 use p007 multipleContfragmentsComplete-template, p007_AAV_03 use p006 multipleContfragmentsComplete-template*

For a quick check of unique fragments 
```bash
sbatch ~/RAAV-30/scripts/UniqueFrag_p007.sh
Rscript ../../scripts/UniqueFrag_p007.R
```
Sort the file by the Reads column
```bash
sort -t, -k1,1 found.p007_AAV_03_frag.csv > sorted_found.p007_AAV_03_frag.csv
# Get unique values based on the Reads column ####
awk -F, '!seen[$1]++' sorted_found.p007_AAV_03_frag.csv > unique_sorted_found.p007_AAV_03_frag.csv
```

## 11.6 RNA count and making .rds files for AAV folders ##
Use p006 multipleContfragmentsComplete.rda
```bash
sbatch ~/RAAV-30/scripts/RNA_count_p007_sub.sh # running in parallel

# Keep individual chunk, applies to p006 (Use p007 multipleContfragmentsComplete.rda)
sbatch ~/RAAV-30/scripts/RNA_count_p006_INDchunk.sh # (took 1 day to finish 30 chunks), use this

# Time consuming part applied, changed to future,future apply, 40 chunks finished in 13 hrs, this is the most efficient for p005
sbatch ~/RAAV-30/scripts/RNA_count_p005_INDchunk_time.sh 
Rscript ../../scripts/RNA_count_p005_INDchunk_time.R
```

## 11.7 RNA count and making .rds files for RAAV-53 and RAAV-72 folders ##
RNA_count (RAAV-53, RAAV-72), p006_tissue use p006 multipleContfragmentsComplete-template, p007_tissue use p007 multipleContfragmentsComplete-template
_p006_
```bash
sbatch ~/RAAV-53/scripts/RNA_count_p006.sh 
Rscript RNA_count_p006_sub.R # (Run using the alignedLibraries.rda, finished in about four hours)
sbatch ~/RAAV-72/scripts/RNA_count_p006.sh 
Rscript RNA_count_p006_sub.R # (Run using the alignedLibraries.rda)

# RNA_count (RAAV-53, RAAV-72) for p005, p007, included alignedLibraries.rda #
sbatch ~/RAAV-72/scripts/RNA_count_p005.sh # use this_run individual sample manually for p005 multipleContfragmentsComplete 
Rscript ../scripts/RNA_count_p005_INDchunk_time.R 

sbatch ~/RAAV-72/scripts/RNA_count_p007_ind_time.sh # use this
Rscript ../scripts/RNA_count_p007_INDchunk_time.R 

sbatch ~/RAAV-53/scripts/RNA_count_p007.sh # use this
Rscript ../scripts/RNA_count_p007_INDchunk_time.R 

sbatch ~/RAAV-53/scripts/RNA_count_p005.sh # use this
Rscript ../scripts/RNA_count_p005_INDchunk_time.R 

sbatch ~/RAAV-72/scripts/RNA_count_p007.sh # Run in bulk (very slow compared to individual manually run)
Rscript RNA_count_p007_sub.R

sbatch ~/RAAV-53/scripts/RNA_count_p007_sub.sh # Run in bulk (very slow compared to individual manually run)
Rscript RNA_count_p007_sub.R
```

## 11.8 Check .rds file and concatenate them properly in each 04_p007_found (RAAV-53 and RAAV-72) ##
```bash
sbatch ~/RAAV-72/scripts/concatenate_foundp007.sh
Rscript scripts/concatenate_found007.R 
```

## 11.9 Collect all .RDS files from different directories and save them into Data_RAAV/02_output folder ##
Run at terminal Data_RAAV 
Use this script for transferring .RDS files from RAAV-30, 60
```bash
Rscript scripts/collectRDS.R 
```

## 11.10 Clean the .RDS file name to be used in normalization ### 
Run at terminal Data_RAAV, apply only for p006
```bash
Rscript scripts/cleanName.R
```

# 12. Data analysis - Normalize Library counts # 
## 12.1 p006 (done in about 1 hour)
```bash
sbatch ~/Data_RAAV/scripts/normalize.sh 
Rscript normalize.R

# To make smaller chunks of p005 completeLibraryRanges 
sbatch ~/Data_RAAV/scripts/chunk.sh # p005 (39 chunks)
Rscript chunk.R

sbatch ~/Data_RAAV/scripts/chunk_p007.sh # p007 (29 chunks)
Rscript chunk_p007.R

# Collect .rds files and make sure not overwrite the existing files with the same name, use this for transferring .RDS files from RAAV-53, 72 
# Remove chunks, p005 ,p006 #
Rscript scripts/collectRDS_existFile.R 
```

## 12.2 To run smaller chunk in a loop of normalization, updates the read count before running the next chunk 
Apply to p005 (RAAV-53; RAAV-72)
```bash
sbatch ~/Data_RAAV/scripts/normalize_42_p005.sh # example for sample 42
Rscript ../scripts/normalize_IND_42_p005.R
sbatch ~/Data_RAAV/scripts/normalize_S59_p005.sh # example for sample S59
Rscript ../scripts/normalize_IND_S59_p005.R
```
## 12.3 To run smaller chunk and save a combined .rds, continue with the step of removing the RNAcount = 1 and normalization 
Apply to p007
```bash
sbatch ~/Data_RAAV/scripts/normalize_42_p007.sh # example for sample 42
Rscript ../scripts/normalize_IND_42_p007.R
sbatch ~/BINP39HooiMin/Data_RAAV/scripts/normalize_S59_p007.sh # example for sample S59
Rscript ../scripts/normalize_IND_S59_p007.R

# This line is important to prevent complain of out of memory: 
Set the future plan with reduced workers
plan(multicore, workers = 8); 
```
## 12.4 Collect normalize .rds files
Collect normalized .rds file and concatenate them into allsamplesDataTable_p007.RDS
```bash
Rscript scripts/collectRDS_normalize.R  

# Collect normalized .rds file and concatenate them into allsamplesDataTable_p005.RDS 
sbatch ~/Data_RAAV/scripts/collectRDS_normalize_concatenatep005.sh
Rscript ../scripts/collectRDS_normalize_concatenatep005.R
```

## 12.5 Pairwise plotting ##
```bash
sbatch ~/Data_RAAV/scripts/PlotGeneCoverage.sh # p006
Rscript ../scripts/PlotAllGenesCoverage_copy2_ind.R # top 25 genes (one plot pairs) # manually pdf done p006

# Use for Pairwise Plotting #
sbatch ~/Data_RAAV/scripts/PlotGeneCoverage_p007_AllVsLib.sh # -- p007 plot mRNA_All vs library top 200 genes (single sample done individually, use RAAV-72 samples, plot72.out)
Rscript ../Data_RAAV/scripts/PlotAllGenesCoverage_copy2_ind_p007_AllVsLib.R

sbatch ~/Data_RAAV/scripts/PlotAccumulate_p007_72.sh # -- p007 Accumulate all RAAV-72 samples, plot mRNA_72 vs library, plot72lib.out
Rscript ../Data_RAAV/scripts/PlotAccumulate_p007_72.R

sbatch ~/Data_RAAV/scripts/PlotAccumulate_p007_72.sh # (accumulate all samples in mRNA_72 vs mRNA_All) - use this
Rscript ../scripts/PlotAccumulate_p007_72.R

sbatch ~/Data_RAAV/scripts/PlotAccumulate_p007_53.sh # (accumulate all samples in mRNA_53 vs mRNA_All) - use this to look for top 100 genes
Rscript ../scripts/PlotAccumulate_p007_53.R
```

## 12.6 Vennplot (Circular barplot) ##
```bash
sbatch ~/Data_RAAV/scripts/VennPlot_p007.sh
Rscript ../scripts/VennPlot_p007.R # (for p007 library)

sbatch ~/Data_RAAV/scripts/PreVenn_p005_53.sh # p005 (done)
Rscript ../scripts/PreVenn_p005_53.R

sbatch ~/Data_RAAV/scripts/PreVenn_p005_72.sh # p005 (need to investigate)
```

## 12.7 Weblogo ##
```bash
sbatch ~/Data_RAAV/scripts/Peptide_cluster.sh # p006
Rscript ../scripts/Peptide_cluster.R

sbatch ~/Data_RAAV/scripts/Peptide_cluster_p007_53.sh # p007 tissue_53
Rscript ../scripts/Peptide_cluster_p007_53.R

sbatch ~/Data_RAAV/scripts/Peptide_cluster_p007_72.sh # p007 tissue_72
Rscript ../scripts/Peptide_cluster_p007_72.R

# Weblogo for selected sample, greedy clustering = 15
sbatch ~/Data_RAAV/scripts/Peptide_cluster_p007_53_SelS.sh # p007 tissue_53
Rscript ../scripts/Peptide_cluster_p007_53_SelS.R
sbatch ~/Data_RAAV/scripts/Peptide_cluster_p007_72_selS.sh # p007 tissue_72
Rscript ../scripts/Peptide_cluster_p007_72_selS.R

# Get top25 motif peptide 
sbatch ~/Data_RAAV/scripts/MostMotif_p007.sh 
Rscript ../scripts/MostMotif_p007.R

# To check if java running correctly, testing at the same directory of the java system_Hammock p007_data 
# Load modules
module load GCCcore/10.3.0
module load GCCcore/11.2.0
module load GCCcore/11.3.0
module load GCCcore/12.3.0
module load ANTLR/2.7.7-Java-11

# Compile the Java code
javac StringLengthChecker.java
javac HHsuiteRunner.java 

# Run the compiled Java class
java StringLengthChecker
java HHsuiteRunner
```



















