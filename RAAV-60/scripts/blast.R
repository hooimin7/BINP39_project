#' ---
#' title: "Blasting CustomArray oligos"
#' author: "Hooi Min Tan Grahn"
#' output: 
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.7in
#' ---

#' This workflow identifies correct fragments from the AAV plasmid library and aligns them to the CustumArray ordered nucleotide fragments using Blastn. Consistant mutations in each fragment/barcode combination are also registered as is the putity of each barcode.

suppressPackageStartupMessages(library(knitr)) 
#+ setup, include=FALSE
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GeneGA))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(biovizBase))
suppressPackageStartupMessages(library(Gviz))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(stringdist))
suppressPackageStartupMessages(library(scales))

opts_chunk$set(fig.width = 7.5, fig.height = 8)
opts_chunk$set(tidy=TRUE)
opts_chunk$set(comment = NA)

strt1<-Sys.time()

#'Load sequences
#'===================
LUT.dna <- fread("fragmentsFinal.csv")

#'Remove constitutive backbone sequences
#'===================
LUT.dna[, Seq := gsub("aacctccagagaggcaacgct","",Seq)]
LUT.dna[, Seq := gsub("gccagacaagcagctaccgca","",Seq)]
LUT.dna[, Seq := toupper(Seq)]
setkey(LUT.dna, "Seq")
LUT.dna <- unique(LUT.dna)
LUT.dna[, Names := Seq]
LUT.dna[, LUTnr := make.names(seq(.N), unique=TRUE)]

# Save the processed data to a CSV file
# fwrite(LUT.dna, file = "LUTdna.csv")

save(LUT.dna,file = "LUTdna.rda")


#' Load the trimmed reads
#' ============================
#+ Loading reads.......

load("LUTdna.rda")

fragments.file <- "../02_pairfq/fragments_paired_reads.fastq.gz"
barcodes.file <- "../02_pairfq/barcodes_paired_reads.fastq.gz"

reads.trim <- readFastq(fragments.file)
reads.BC <- readFastq(barcodes.file)

#' Make CustomArray reference index for Blast
#' ============================
#+ Making Bowtie index.......

LUT.fa <- tempfile(pattern = "LUT_", tmpdir = tempdir(), fileext = ".fa")
LUT.seq = ShortRead(DNAStringSet(LUT.dna$Seq), BStringSet(1:length(LUT.dna$LUTnr)))
writeFasta(LUT.seq,LUT.fa)


print(LUT.fa)

# Read the file
# fasta_content <- readLines(LUT.fa)

# Print the content
# print(fasta_content)

#'Save unique fragments as fasta file
#'===================
unique.reads <- unique(sread(reads.trim))

#'Select subset
#'===================
# Define the directory path
# dir_path <- "/lunarc/nobackup/projects/lu2024-17-19/RAAV-60/p006/03_blast"

# unique.reads <- ShortRead(DNAStringSet(unique.reads), BStringSet(1:length(unique.reads)))

# # Use the directory path in the tempfile() function
# fragments_unique_fa <- tempfile(pattern = "FragUnique_", tmpdir = dir_path, fileext = ".fa")
# writeFasta(unique.reads,fragments.unique.fa)

unique.reads <- ShortRead(DNAStringSet(unique.reads), BStringSet(1:length(unique.reads)))
fragments.unique.fa <- tempfile(pattern = "FragUnique_", tmpdir = tempdir(), fileext = ".fa")
writeFasta(unique.reads,fragments.unique.fa)

# Print the path to the temporary file
print(fragments.unique.fa)

# Read the file
# fasta_content_sub <- readLines(fragments.unique.fa)

# Print the content
# print(fasta_content_sub)

#'Align against the library using blast
#'===================

blast.db <- tempfile(pattern = "blastDB_", tmpdir = tempdir(), fileext = ".db")
blast.out <- tempfile(pattern = "blastOut_", tmpdir = tempdir(), fileext = ".txt")

sys.out <-  system(paste("makeblastdb -in ", LUT.fa,
                         " -out ",blast.db," -dbtype nucl -title LUT -parse_seqids 2>&1",  sep = ""),
                   intern = TRUE, ignore.stdout = FALSE) 

sys.out <- as.data.frame(sys.out)

print(sys.out)

colnames(sys.out) <- c("blastn database generation")
invisible(sys.out[" "] <- " ")
knitr::kable(sys.out[1:(nrow(sys.out)),], format =  "latex", booktabs = T) %>%
  kable_styling(latex_options = "striped")
#   as.character() %>%
#   htmltools::HTML() %>%
#   htmltools::browsable()

sys.out <-  system(paste("export SHELL=/bin/sh; cat ",fragments.unique.fa," | parallel --block ",
                         floor(length(unique.reads)/detectCores()),
                         " --recstart '>' --pipe 'blastn -max_target_seqs 25 -word_size 7",
                         " -num_threads 1 -outfmt 10 -db ", blast.db,
                         " -query - '> ", blast.out, " 2>&1",  sep = ""),
                   intern = TRUE, ignore.stdout = FALSE) 


system(paste("gzip -c ", blast.out, " > ./data/blastOutput.csv.gz", sep=" "))

table.blastn <- data.table(scan(file="./data/blastOutput.csv.gz", what="character", sep=";") , keep.rownames=FALSE, key="V1")

names(table.blastn)

if (length(grep("Warning",table.blastn$V1)) != 0) {
  warnings.out <- unique(table.blastn[grep("Warning",table.blastn$V1),])
  table.blastn <- table.blastn[-grep("Warning",table.blastn$V1),]
  setnames(warnings.out,"V1", c("blastn Warnings"))
  # knitr::kable(warnings.out[1:(nrow(warnings.out)),], format = "markdown")
}

table.blastn[,c("Reads","Sequence","identity","alignmentLength","mismatches",
                "gapOpens", "q_start", "q_end", "s_start", "s_end", 
                "evalue","bitScore") := tstrsplit(V1,",",fixed=TRUE),]


table.blastn[,Reads:= as.character(sread(unique.reads[as.integer(Reads)]))]
print(names(table.blastn))
table.blastn[,Sequence:= as.character(sread(LUT.seq[as.integer(Sequence)]))]
print(LUT.seq)
setkey(table.blastn,Sequence)
setkey(LUT.dna,Seq)
print(names(LUT.dna))
table.blastn<- table.blastn[LUT.dna, nomatch=0]
table.blastn[,c("V1","identity","alignmentLength","gapOpens", "q_start", 
                "q_end", "s_start", "s_end", "evalue","Sequence","Names"):=NULL]
gc() #garbage collection to reduce memory foot print. Can be removed for speed

table.blastn[,bitScore:= as.numeric(bitScore)]
table.blastn[,mismatches:= as.numeric(mismatches)]

setkeyv(table.blastn,c("Reads","LUTnr"))
setorder(table.blastn,Reads,LUTnr,-bitScore) #This makes sure that a fragment is only aligned once to the reference in the top ten matches
table.blastn <- unique(table.blastn, by=c("Reads","LUTnr"))

gc() #garbage collection to reduce memory foot print. Can be removed for speed

table.blastn.topHit <- table.blastn[table.blastn[, .I[which.max(bitScore)], by="Reads"]$V1] # Select only rows with the highest bitScore

full.table <- data.table(Reads=as.character(sread(reads.trim)),
                         BC=as.character(sread(reads.BC)),
                         key="Reads")
all.reads <- nrow(full.table)

full.table <- full.table[table.blastn.topHit,nomatch=0] # Merge reads with the top hit alignment

print(paste("Alignment percentage:", percent(nrow(full.table)/all.reads)))
