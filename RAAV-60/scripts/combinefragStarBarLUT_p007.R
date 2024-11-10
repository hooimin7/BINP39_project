rm(list=ls()) # Clear workspace
library(data.table)
library(utils)
library(ShortRead)
library(digest)
library(scales)
library(Biostrings)
library(dplyr)

# Load sequences
# ===================
# Read LUT.dna CSV file
LUT.dna <- fread("LUTdna.csv")

# Create ShortRead object
LUT.seq <- ShortRead(DNAStringSet(LUT.dna$Seq), BStringSet(1:length(LUT.dna$LUTnr)))

# File paths for fragments and barcodes
fragments.file <- "../06_library/fragments_7_FragBar_paired_again.fastq.gz"
barcodes.file <- "../06_library/barcodes_7_FragBar_paired_again.fastq.gz"

# Read FASTQ files
reads.trim <- readFastq(fragments.file)
reads.BC <- readFastq(barcodes.file)

# Extract sequences from reads.trim
sequences <- as.character(sread(reads.trim))

# Count the occurrences of each sequence
sequence_counts <- table(sequences)

# Convert to a data.table
sequence_counts.dt <- as.data.table(sequence_counts, keep.rownames = "Seq")

# Rename columns for clarity
setnames(sequence_counts.dt, c("Seq", "Frequency"))

# Convert LUT.dna to data.table if not already
setDT(LUT.dna)

# Set keys
setkey(LUT.dna, "Seq")
setkey(sequence_counts.dt, "Seq")

# Perform the join
frag.LUT <- LUT.dna[sequence_counts.dt, nomatch = 0]


print("frag.LUT join")
print(frag.LUT)

# Set keys and order for unique alignment
setkeyv(frag.LUT, c("Seq", "LUTnr"))
setorder(frag.LUT, Seq, LUTnr, -Frequency)

# Ensure unique alignment
frag.LUT <- unique(frag.LUT, by = c("Seq", "LUTnr"))

# Create full table with reads and barcodes
full.table <- data.table(Reads = as.character(sread(reads.trim)), BC = as.character(sread(reads.BC)),
                         key = "Reads")

# Reduce barcode
out.name.BC.star <- "../07_starcode/01_data/barcodes_7_reduced.txt"

table.BC.sc <- data.table(read.table(out.name.BC.star, header = FALSE, row.names = 1,
                                     skip = 0, sep = "\t", stringsAsFactors = FALSE, fill = FALSE), keep.rownames = TRUE,
                          key = "rn") #, nrows = 1000

table.BC.sc[, `:=`(V2, NULL)]

table.BC.sc <- table.BC.sc[, strsplit(as.character(V3), ",", fixed = TRUE),
                           by = rn]

unique_reads_BC <- length(unique(sread(reads.BC)))
unique_table_BC_sc <- length(unique(table.BC.sc$V1) %in% unique(sread(reads.BC)))
SC.droppedBC <- unique_reads_BC - unique_table_BC_sc

setnames(table.BC.sc, c("V1", "rn"), c("BC", "scBC"))

# Replacing barcodes with Starcode reduced versions
# Setting keys for data tables
setkey(full.table, BC)
setkey(table.BC.sc, BC)

# Merging data tables
full.table <- full.table[table.BC.sc, nomatch = 0]

print("full.table")
print(head(full.table))

# Renaming columns
setnames(full.table, c("BC", "scBC"), c("oldBC", "BC"))
setkey(full.table, BC)

# Calculating unique barcodes
RetainedBC <- length(unique(full.table$oldBC))
scBC <- length(unique(full.table$BC))
print(paste("Original unique barcodes:", RetainedBC))
print(paste("SC reduced unique barcodes:", scBC))

# Creating and formating tables
table.frag.oldBC <- data.table(as.data.frame((rev(sort(table(full.table$oldBC))))[1:10],
                                             row.names = "Var1"), keep.rownames = TRUE)

# In R versions below 3.3 remove, row.names = 'Var1' to make this compatible

setnames(table.frag.oldBC, colnames(table.frag.oldBC), c("Original BC", "Count"))

table.frag <- data.table(as.data.frame((rev(sort(table(full.table$BC))))[1:10],
                                       row.names = "Var1"), keep.rownames = TRUE)
# In R versions below 3.3 remove, row.names = 'Var1' to make this compatible
setnames(table.frag, colnames(table.frag), c("SC reduced BC", "Count"))

# Removing oldBC column
invisible(full.table[, `:=`(oldBC, NULL)])

# Ordering the table
full.table <- full.table[order(full.table$BC), ]

# Splitting reads into single-read and multi-read barcodes
temp.table.single <- full.table[full.table[, .I[.N == 1], by = "BC"]$V1] # BC appears once
temp.table.multi <- full.table[full.table[, .I[.N > 1], by = "BC"]$V1] # BC appears more than once

# Adding columns to single read table
temp.table.single[, `:=`(c("mCount", "tCount"), 1)] # Add mCount and tCount col with a vaule of 1
temp.table.single$Mode <- "Amb" # Add a Mode column with value (Ambiguos)

print("temp.table.single")
print(temp.table.single)

# Setting key for multireads table
setkeyv(temp.table.multi, c("BC"))

# Aggregating multi-read table
temp.table.multi[, `:=`(tCount, .N), by = key(temp.table.multi)] # Set tCount to the number of rows in each group

temp.table.multi$Mode <- "Def" # Add definitive column
temp.table.multi <- unique(temp.table.multi) # Remove duplicate rows

print("temp.table.multi")
print(temp.table.multi)

print("Utilized reads.......")
print(nrow(full.table))
print("Whereof single reads.......")
print(nrow(temp.table.single))

# Splitting multi-read barcodes into clean and chimeric
setkeyv(temp.table.multi, "BC")

# Filter for barcodes
temp.table.multi.clean <- temp.table.multi[temp.table.multi[, .I[.N == 1], by = "BC"]$V1] # BC appears only once

# Filter for non-unique barcodes
temp.table.multi <- temp.table.multi[temp.table.multi[, .I[.N > 1], by = "BC"]$V1] # BC appears more than once

# Update mCount column
temp.table.multi.clean[, `:=`(mCount, tCount)]

print("temp.table.multi.clean")
print(temp.table.multi.clean)

print("Clean multi-read barcodes.......")
print(nrow(temp.table.multi.clean))
print("Chimeric multi-read barcodes.......")
print(length(unique(temp.table.multi$BC)))

# Calculate consensus alignment of chimeric barcodes
setkey(temp.table.multi, "BC")

# Update mCount column
temp.table.multi[, mCount := tCount]

# Aggregate tCount by BC
temp.table.multi[, tCount := sum(tCount), by = "BC"]

setkey(temp.table.multi, "Reads")

# Change name in data.table
setnames(frag.LUT, old = "Seq", new = "Reads")

setkey(frag.LUT, "Reads")

# Remove names from frag.LUT
frag.LUT[, `:=`(Names, NULL)]

# Merge temp.table.multi with table.blastn using by = .EACHI to avoid large allocation
temp.table.multi <- temp.table.multi[frag.LUT, nomatch = 0, allow.cartesian = TRUE]

# Remove Type from temp.table.multi
temp.table.multi[, `:=`(Type, NULL)]

print("after merge temp.table.multi and table.blastn")
print(temp.table.multi)

setkeyv(temp.table.multi, c("BC", "LUTnr"))

# Aggregate column by key
temp.table.multi[, `:=`(mCount, sum(mCount)), by = key(temp.table.multi)]

# Remove duplicates from temp.table.multi
temp.table.multi <- unique(temp.table.multi, by = c("BC", "LUTnr"))

setkeyv(temp.table.multi, "BC")

# Select rows with the highest mCount
temp.table.multi <- temp.table.multi[temp.table.multi[, .I[mCount == max(mCount)],
                                                      by = key(temp.table.multi)]$V1]

# Set Mode to "Amb" for mCount = 1
temp.table.multi[temp.table.multi$mCount == 1, Mode := "Amb"]

print(paste("Number of barcodes with false mCount:", nrow(temp.table.multi[mCount > tCount])))

# Combine temp.table.multi and temp.table.multi.clean
temp.table.multi.consensus <- bind_rows(temp.table.multi, temp.table.multi.clean)

print(paste("Total number of definitive Barcodes:", length(grep("Def", temp.table.multi.consensus$Mode))))
print(paste("Total number of ambiguous Barcodes:", length(grep("Amb", temp.table.multi.consensus$Mode))))
print(paste("Total number of single-read Barcodes:", nrow(temp.table.single)))

print("temp.table.multi.consensus")
print(temp.table.multi.consensus)
print(paste("Number of rows in temp.table.consensus:", nrow(temp.table.multi.consensus)))

print("temp.table.multi at the end")
print(temp.table.multi)
print(paste("Number of rows in temp.table.multi:", nrow(temp.table.multi)))

print("temp.table.multi.clean at the end")
print(temp.table.multi.clean)
print(paste("Number of rows in temp.table.multi.clean:", nrow(temp.table.multi.clean)))

# Check for duplicates in temp.table.multi.consensus
num_duplicates <- nrow(temp.table.multi.consensus) - nrow(unique(temp.table.multi.consensus))
print(paste("Number of duplicate rows in temp.table.multi.consensus:", num_duplicates))

# Verify the number of unique definitive barcodes
unique_def_barcodes <- length(unique(temp.table.multi.consensus[Mode == "Def", BC]))
print(paste("Unique definitive barcodes:", unique_def_barcodes))

# Combine all data tables
output.Table <- bind_rows(temp.table.multi.consensus, temp.table.single)

# Save the output table as CSV
write.csv(output.Table, file = "02_analysis/multipleContfragmentsComplete.csv", row.names = FALSE)

# Save the output table as RDA
save(output.Table, file = "02_analysis/multipleContfragmentsComplete.rda")


##################################

