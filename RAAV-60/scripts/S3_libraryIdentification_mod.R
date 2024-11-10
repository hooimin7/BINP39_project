
#' ---
#' title: "Reverse mapping of CustumArray oligos to original proteins"
#' author: "Tomas Bjorklund"
#' edited by: "Jaro Steindorff and Hooi Min"
#' output: alignedLibraries.rda -- a table of the fragments with their origin protein and the LUTnr, strucutre and sequence.
#'         This workflow identifies correct fragments from the Cre-recombined AAV plasmid library and aligns them to the CustumArray ordered nucleotide fragments using Blastn.
#'         Consistant mutations in each fragment/barcode combination are also registered as is the putity of each barcode.
#' ---

#  This will make R-functions such as library() and install.packages() use this directory:
.libPaths(c('~/MyRextensions', .libPaths()))

#+ setup, include=FALSE
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(GeneGA))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(Biostrings))

source(file.path("../../scripts/GeneCodon.R"))
#Override the GeneCodon function with local version containing human codons
unlockBinding("GeneCodon", as.environment("package:GeneGA"))
assign("GeneCodon", GeneCodon, as.environment("package:GeneGA"))

#'Load sequences
#'===================
LUT.dna <- fread("../04_blast/LUTdna.csv")
LUT.dna <- data.table(LUT.dna)

cat("First few rows of LUT.dna:\n")
print(head(LUT.dna))

# Extract sequences and IDs
lut_sequences <- LUT.dna$Seq
ids <- LUT.dna$LUTnr

# Debug: Print the extracted sequences and IDs
cat("Extracted sequences:\n")
print(head(lut_sequences))
cat("Extracted IDs:\n")
print(head(ids))

#' Save fasta files for Bowtie alignments
#' ===================

# Create a ShortRead object for LUT.7aa sequences
LUT.7.seq <- ShortRead(DNAStringSet(lut_sequences), BStringSet(ids))

cat("ShortRead object LUT.7.seq:\n")
print(LUT.7.seq)

# Write the LUT.7aa sequences to a fasta file
LUT.7aa.fa <- "LUT_7aa.fa"
writeFasta(LUT.7.seq,LUT.7aa.fa)

#'Build Bowtie index
#'===================
# Read the TSV file
tsv_file <- "20220904_Sequences_for_fragmentation_uppercase.tsv"
seqs.TSV <- read.delim(tsv_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Extract gene names and sequences
gene_names <- seqs.TSV$Genename
sequences <- seqs.TSV$Sequence

# Create a DNAStringSet object
dna_sequences <- DNAStringSet(sequences)
names(dna_sequences) <- gene_names

# Write to FASTA file
fasta_file <- "Sequences_for_fragmentation.fa"
writeXStringSet(dna_sequences, filepath = fasta_file)

output_file <- 'Sequences_for_fragmentation_mod.fa'
cleaned_output_file <- 'Sequences_for_fragmentation_cleaned.fa'

# Read the input file
lines <- readLines(fasta_file)

# Initialize a counter for numbering the sequences
counter <- 1

# Open the output file for writing
outfile <- file(output_file, 'w')

# Process each line
for (line in lines) {
  if (startsWith(line, '>')) {
    header <- trimws(line)
    new_header <- paste0(header, ",No,,#", counter, substring(header, 2))
    writeLines(new_header, outfile)
    counter <- counter + 1
  } else {
    writeLines(line, outfile)
  }
}

# # Process each line
# for (line in lines) {
#   if (startsWith(line, '>')) {
#     header <- trimws(line)
#     new_header <- paste0(header, ",No,,#", counter, ",")
#     writeLines(new_header, outfile)
#     counter <- counter + 1
#   } else {
#     writeLines(line, outfile)
#   }
# }

# Close the output file
close(outfile)

# read in the sequneces from the original fasta file
seqs.original <- readDNAStringSet(output_file)

# Extract IDs from the headers (e.g., #1, #2, etc.)
seq_ids <- sapply(names(seqs.original), function(header) {
  strsplit(header, ",")[[1]][4]
})

# Function to clean sequences by replacing invalid characters
clean_sequence <- function(seq) {
  seq <- as.character(seq)
  seq <- gsub("N", "A", seq)  # Replace 'N' with 'A'
  seq <- gsub("[^ATCG]", "A", seq)  # Replace any non-ATCG characters with 'A'
  return(DNAString(seq))
}


# Apply the cleaning function to all sequences
cleaned_sequences <- DNAStringSet(sapply(seqs.original, clean_sequence))

# Set the names of the cleaned sequences to the extracted IDs
names(cleaned_sequences) <- seq_ids

# Write the cleaned sequences back to a new FASTA file
writeXStringSet(cleaned_sequences, filepath = cleaned_output_file)

print("Cleaned_seq")
print(cleaned_sequences)

# Translate the cleaned sequences to amino acids
seqs.AA <- Biostrings::translate(cleaned_sequences, genetic.code=GENETIC_CODE, if.fuzzy.codon="error")

print("Translated aa seq")
print(seqs.AA)

# Handle stop codons in the translated sequences
handle_stop_codons <- function(aa_seq) {
  aa_seq <- as.character(aa_seq)
  aa_seq <- gsub("\\*", "", aa_seq)  # Remove '*' (stop codons)
  return(AAString(aa_seq))
}

# Apply the function to handle stop codons
handled_seqs.AA <- AAStringSet(sapply(seqs.AA, handle_stop_codons))

# Print the handled amino acid sequences
print("Handled amino acid sequences:")
print(handled_seqs.AA)

# load the function to convert amino acids to DNA
source("../../scripts/scAAtoDNA.R")

# Convert amino acids to DNA
optimized_dna_sequences <- DNAStringSet(sapply(handled_seqs.AA, function(x) AAtoDNA(x, species="hsa")))

# Create a BStringSet object for the IDs
seq_ids_bstring <- BStringSet(gsub("([ ])", "_", seq_ids))

# Create a ShortRead object with the optimized sequences and IDs
seqs.optimized <- ShortRead(optimized_dna_sequences, seq_ids_bstring)

# write the optimized sequences to a fasta file
bowtie.fasta <- "bowtie.fa"
writeFasta(seqs.optimized,bowtie.fasta)

# create a temporary file for the bowtie index
# bowtie.idx <- tempfile(pattern = "IDX_bowtie_", tmpdir = tempdir(), fileext = "")
bowtie.idx <- tempfile(pattern = "IDX_bowtie_", tmpdir = ".", fileext = "")

sys.out <-  system(paste("bowtie2-build",bowtie.fasta,bowtie.idx, "2>&1",  sep = " "), 
                   intern = TRUE, ignore.stdout = FALSE) 

# ' Align fragments to reference
# ' ============================

# Function to align sequences and process the results
align_and_process <- function(fa_file, lut, bowtie_idx) {
  # Create a temporary file name for Bowtie output
  name.bowtie <- tempfile(pattern = "bowtie_", tmpdir = tempdir(), fileext = "")

  # Run Bowtie2 alignment
  sys.out <- system(paste("bowtie2 --non-deterministic --threads ", detectCores(),
                          " --very-sensitive -f -a",
                          " -x ", bowtie_idx, " -U ", fa_file, " -S ", 
                          name.bowtie, ".sam 2>&1", sep = ""), 
                    intern = TRUE, ignore.stdout = FALSE)

                    # Run samtools to convert the SAM file to a BAM file
  command.args <- paste("view -@ ", detectCores(), " -Sb ", name.bowtie, ".sam > ",
                        name.bowtie, ".bam", sep = "")
  system2("samtools", args = command.args, stdout = TRUE, stderr = TRUE)
  
  # Run samtools to sort the BAM file
  command.args <- paste("sort -@ ", detectCores(), " ", name.bowtie, ".bam -o ",
                        name.bowtie, "_sort.bam", sep = "")
  system2("samtools", args = command.args, stdout = TRUE, stderr = TRUE)

  # Read the sorted BAM file into a GAlignments object
  frag_ranges <- readGAlignments(paste(name.bowtie, "_sort.bam", sep = ""), use.names = TRUE)
  
  # Return the results as a list
  list(
    frag_ranges = frag_ranges,
    total = length(names(frag_ranges)),
    unique = length(unique(names(frag_ranges))),
    lut_unique = length(unique(lut$Sequence))
  )
}

results_7aa <- align_and_process(LUT.7aa.fa, LUT.dna, bowtie.idx)

print(head(results_7aa))

# Merge all the aligned sequences into one object
allFragments.ranges <- results_7aa$frag_ranges

print(head(allFragments.ranges))

# Annotate the merged sequences with their LUT number
mcols(allFragments.ranges)$LUTnr <- names(allFragments.ranges)
# Set the key for the LUT data frame
setkey(LUT.dna, LUTnr)
# Annotate the merged sequences with their corresponding sequence from the LUT
mcols(allFragments.ranges)$Sequence <- LUT.dna[mcols(allFragments.ranges)$LUTnr]$Seq

# Save the merged and annotated sequences to a file
save(allFragments.ranges, file = "alignedLibraries.rda")

# Save the merged and annotated sequences to a csv file
write.csv(allFragments.ranges, file = "alignedLibraries.csv", row.names = FALSE)


devtools::session_info()

