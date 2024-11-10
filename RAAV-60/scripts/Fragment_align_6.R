rm(list=ls()) # Clear workspace

library(dplyr)
library(data.table)
library(utils)
library(ShortRead)
library(digest)
library(scales)
library(Biostrings)

# Load sequences
# ===================
# Read LUT.dna CSV file
LUT.dna <- fread("LUTdna.csv")

# Create ShortRead object
LUT.seq <- ShortRead(DNAStringSet(LUT.dna$Seq), BStringSet(1:length(LUT.dna$LUTnr)))

# File paths for fragments and barcodes
fragments.file <- "../06_library/fragments_6_FragBar_paired_again.fastq.gz"
barcodes.file <- "../06_library/barcodes_6_FragBar_paired_again.fastq.gz"

# Read FASTQ files
reads.trim <- readFastq(fragments.file)
reads.BC <- readFastq(barcodes.file)

# Extract the ID and sequence lines
ids <- as.character(id(reads.trim))
sequences <- as.character(sread(reads.trim))

# Create a data frame with the ID in the first column and the sequence in the second column
fragment_data <- data.frame(ID = ids, Sequence = sequences, stringsAsFactors = FALSE)

# Count the occurrences of each sequence
fragment_counts <- table(sequences)

# Convert to a data.table
sequence_counts.dt <- as.data.table(fragment_counts, keep.rownames = "Sequence")

# Rename columns for clarity
setnames(sequence_counts.dt, c("Sequence", "Frequency"))

# Merge the count information with the main fragment data
fragment_data <- merge(fragment_data, sequence_counts.dt, by = "Sequence", all.x = TRUE)

# Extract the relevant part of the ID in fragment_data
fragment_data$ID <- sub(" .*", "", fragment_data$ID)

# Read blast output into data.table
table.blastn <- data.table(scan(file = "01_blastout/blast_output_6.csv.gz", what = "character",
                                sep = ";"), keep.rownames = FALSE, key = "V1")

# Handle warnings in blast output
if (length(grep("Warning", table.blastn$V1)) != 0) {
  warnings.out <- unique(table.blastn[grep("Warning", table.blastn$V1), ])
  table.blastn <- table.blastn[-grep("Warning", table.blastn$V1), ]
  setnames(warnings.out, "V1", c("blastn Warnings"))
}

# Split V1 column into multiple columns
table.blastn[, `:=`(c("Reads", "Sequence", "identity", "alignmentLength", "mismatches",
                      "gapOpens", "q_start", "q_end", "s_start", "s_end", "evalue", "bitScore"),
                    tstrsplit(V1, ",", fixed = TRUE)), ]

# Remove V1 column
table.blastn[, V1 := NULL]

# Rename the "Reads" column to "Id" and the "Sequence" column to "LUTnr"
setnames(table.blastn, old = c("Reads", "Sequence"), new = c("ID", "LUTnr"))

# Convert data frames to data.table objects
setDT(table.blastn)
setDT(fragment_data)

# Set keys for merging
setkey(table.blastn, ID)
setkey(fragment_data, ID)

# Merge table.blastn with fragment_data 
table.blastn <- table.blastn[fragment_data, nomatch = 0, allow.cartesian = TRUE]

# Remove ID column
table.blastn[, ID := NULL]

# Set keys for merging
setkey(table.blastn, LUTnr)
setkey(LUT.dna, LUTnr)

# Merge table.blastn with LUT.dna
table.blastn <- table.blastn[LUT.dna, nomatch = 0, allow.cartesian = TRUE]

# Remove unnecessary columns
table.blastn[, `:=`(c("identity", "alignmentLength", "gapOpens", "q_start",
                      "q_end", "s_start", "s_end", "evalue"), NULL)]

# Convert bitScore and mismatches to numeric
table.blastn[, `:=`(bitScore, as.numeric(bitScore))]
table.blastn[, `:=`(mismatches, as.numeric(mismatches))]

# Set keys and order for unique alignment
setkeyv(table.blastn, c("Sequence", "LUTnr"))
setorder(table.blastn, Sequence, LUTnr, -bitScore) 

# Ensure unique alignment
table.blastn <- unique(table.blastn, by = c("Sequence", "LUTnr"))

# Get top hit alignment
table.blastn.topHit <- table.blastn[table.blastn[, .I[which.max(bitScore)],
                                                 by = "Sequence"]$V1]


# Create full table with reads and barcodes
full.table <- data.table(Sequence = as.character(sread(reads.trim)), BC = as.character(sread(reads.BC)),
                         key = "Sequence")

# Get total number of reads
all.reads <- nrow(full.table)

# Merge full.table with top hit alignment
full.table <- full.table[table.blastn.topHit, nomatch = 0] # Merge reads with the top hit alignment

# print("full.table after tophit")
# print(full.table)

# Print alignment percentage
print(paste("Alignment percentage:", percent(nrow(full.table)/all.reads)))

# Save the output table as CSV
write.csv(full.table, file = "02_analysis/fragment_align_6.csv", row.names = FALSE)

# Read the CSV file without the header
fragment_align <- fread("02_analysis/fragment_align_6.csv", header = TRUE)

# Extract the header
header <- colnames(fragment_align)

# Sort the data by the bitScore and Frequency columns in descending order, skipping the header
sorted_fragment_align <- fragment_align %>%
  arrange(desc(bitScore), desc(Frequency))

# Add the header back
colnames(sorted_fragment_align) <- header

# Save the sorted data to a new CSV file
write.csv(sorted_fragment_align, file = "02_analysis/sorted_fragment_align_6.csv", row.names = FALSE)

DNA_pscAAVlib <- full.table

save(DNA_pscAAVlib, file = "02_analysis/scBC_DNA_pscAAVlib_p006.rda")


# Zipping the output table
file_to_zip <- "02_analysis/sorted_fragment_align_6.csv"

# Define the name of the zip file
zip_file <- "02_analysis/sorted_fragment_align_6.zip"

# Zip the file
zip(zip_file, file_to_zip)

# Save the output table as RDA
save(full.table, file = "02_analysis/sorted_fragment_align_6.rda")



