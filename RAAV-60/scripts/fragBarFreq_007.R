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

# Assuming temp.table.single is already defined and contains the single-read barcodes

# Count the number of reads for each clean barcode
single_read_barcodes <- temp.table.single[, .N, by = BC]
setnames(single_read_barcodes, "N", "BcCount")

# Merge the counts back with the original data to trace the reads
single_traced_reads <- merge(temp.table.single, single_read_barcodes, by = "BC")

# Count the occurrences of each read
single_read_counts <- single_traced_reads[, .N, by = Reads]
setnames(single_read_counts, "N", "ReadCount")

# Merge the read counts back with the traced reads data
single_traced_reads <- merge(single_traced_reads, single_read_counts, by = "Reads")

# Classify reads as unique or repetitive
single_traced_reads[, ReadType := ifelse(ReadCount == 1, "Unique", "Repetitive")]

# Summarize the total number of unique and repetitive reads
summary_counts_single <- single_traced_reads[, .N, by = ReadType]
setnames(summary_counts_single, "N", "Count")

print("Summary counts single") 
print(summary_counts_single)

# Filter out repetitive reads into another variable
single_repetitive_reads <- single_traced_reads[ReadType == "Repetitive"]

# Keep unique reads in the same variable
single_traced_reads <- single_traced_reads[ReadType == "Unique"]

# Remove BcCount, ReadCount, ReadType from single_traced_reads
single_traced_reads[, `:=`(c("BcCount", "ReadCount", "ReadType"), NULL)]

# # Save the unique and repetitive reads to CSV files
# write.csv(single_traced_reads, file = "unique_reads.csv", row.names = FALSE)
# write.csv(single_repetitive_reads, file = "repetitive_reads.csv", row.names = FALSE)
# 
# # Save the unique and repetitive reads to RDA files
# save(single_traced_reads, file = "single_unique_reads.rda")
# save(single_repetitive_reads, file = "single_repetitive_reads.rda")

# Setting key for multireads table
setkeyv(temp.table.multi, c("BC"))

# Aggregating multi-read table
temp.table.multi[, `:=`(tCount, .N), by = key(temp.table.multi)] # Set tCount to the number of rows in each group

temp.table.multi$Mode <- "Def" # Add definitive column
temp.table.multi <- unique(temp.table.multi) # Remove duplicate rows

print("Utilized reads.......")
print(nrow(full.table))
print("Whereof single reads.......")
print(nrow(single_traced_reads))

# Splitting multi-read barcodes into clean and chimeric
setkeyv(temp.table.multi, "BC")

# Filter for barcodes
temp.table.multi.clean <- temp.table.multi[temp.table.multi[, .I[.N == 1], by = "BC"]$V1] # BC appears only once

# Assuming temp.table.multi.clean is already defined and contains the clean multi-read barcodes

# Count the number of reads for each clean barcode
multi_read_clean_barcodes <- temp.table.multi.clean[, .N, by = BC]
setnames(multi_read_clean_barcodes, "N", "BcCount")

# Merge the counts back with the original data to trace the reads
traced_reads <- merge(temp.table.multi.clean, multi_read_clean_barcodes, by = "BC")

# Count the occurrences of each read
read_counts <- traced_reads[, .N, by = Reads]
setnames(read_counts, "N", "ReadCount")

# Merge the read counts back with the traced reads data
traced_reads <- merge(traced_reads, read_counts, by = "Reads")

# Classify reads as unique or repetitive
traced_reads[, ReadType := ifelse(ReadCount == 1, "Unique", "Repetitive")]

# Summarize the total number of unique and repetitive reads
summary_counts <- traced_reads[, .N, by = ReadType]
setnames(summary_counts, "N", "Count")

print("Summary counts clean") 
print(summary_counts)

# Filter out repetitive reads into another variable
repetitive_reads <- traced_reads[ReadType == "Repetitive"]

# Keep unique reads in the same variable
traced_reads <- traced_reads[ReadType == "Unique"]

# Remove BcCount, ReadCount, ReadType from traced_reads
traced_reads[, `:=`(c("BcCount", "ReadCount", "ReadType"), NULL)]

# # Save the unique and repetitive reads to CSV files
# write.csv(traced_reads, file = "unique_reads.csv", row.names = FALSE)
# write.csv(repetitive_reads, file = "repetitive_reads.csv", row.names = FALSE)
# 
# # Save the unique and repetitive reads to RDA files
# save(traced_reads, file = "unique_reads.rda")
# save(repetitive_reads, file = "repetitive_reads.rda")

# Filter for non-unique barcodes
temp.table.multi <- temp.table.multi[temp.table.multi[, .I[.N > 1], by = "BC"]$V1] # BC appears more than once

# Assuming temp.table.multi n is already defined and contains the multi-read barcodes

# Count the number of reads for each clean barcode
multi_read_barcodes <- temp.table.multi[, .N, by = BC]
setnames(multi_read_barcodes, "N", "BcCount")

# Merge the counts back with the original data to trace the reads
multi_traced_reads <- merge(temp.table.multi, multi_read_barcodes, by = "BC")

# Count the occurrences of each read
multi_read_counts <- multi_traced_reads[, .N, by = Reads]
setnames(multi_read_counts, "N", "ReadCount")

# Merge the read counts back with the traced reads data
multi_traced_reads <- merge(multi_traced_reads, multi_read_counts, by = "Reads")

# Classify reads as unique or repetitive
multi_traced_reads[, ReadType := ifelse(ReadCount == 1, "Unique", "Repetitive")]

# Summarize the total number of unique and repetitive reads
summary_counts_multi <- multi_traced_reads[, .N, by = ReadType]
setnames(summary_counts_multi, "N", "Count")

print("Summary counts multi") 
print(summary_counts_multi)

# Filter out repetitive reads into another variable
multi_repetitive_reads <- multi_traced_reads[ReadType == "Repetitive"]

# Keep unique reads in the same variable
multi_traced_reads <- multi_traced_reads[ReadType == "Unique"]

# Remove BcCount, ReadCount, ReadType from traced_reads
multi_traced_reads[, `:=`(c("BcCount", "ReadCount", "ReadType"), NULL)]

# # Save the unique and repetitive reads to CSV files
# write.csv(multi_traced_reads, file = "unique_reads.csv", row.names = FALSE)
# write.csv(multi_repetitive_reads, file = "repetitive_reads.csv", row.names = FALSE)
# 
# # Save the unique and repetitive reads to RDA files
# save(multi_traced_reads, file = "multi_unique_reads.rda")
# save(multi_repetitive_reads, file = "multi_repetitive_reads.rda")

# # Update mCount column
# temp.table.multi.clean[, `:=`(mCount, tCount)]

# Update mCount column
traced_reads[, `:=`(mCount, tCount)]

# print("Clean multi-read barcodes.......")
# print(nrow(temp.table.multi.clean))

print("Clean multi-read barcodes.......")
print(nrow(traced_reads))

# print("Chimeric multi-read barcodes.......")
# print(length(unique(temp.table.multi$BC)))

print("Chimeric multi-read barcodes.......")
print(length(unique(multi_traced_reads$BC)))

# # Calculate consensus alignment of chimeric barcodes
# setkey(temp.table.multi, "BC")

# Calculate consensus alignment of chimeric barcodes
setkey(multi_traced_reads, "BC")

print("Multi-read barcodes.......")
print(nrow(multi_traced_reads))

# # Update mCount column
# temp.table.multi[, mCount := tCount]

# Update mCount column
multi_traced_reads[, mCount := tCount]

# # Aggregate tCount by BC
# temp.table.multi[, tCount := sum(tCount), by = "BC"]

# Aggregate tCount by BC
multi_traced_reads[, tCount := sum(tCount), by = "BC"]

# setkey(temp.table.multi, "Reads")

setkey(multi_traced_reads, "Reads")

# Change name in data.table
setnames(frag.LUT, old = "Seq", new = "Reads")

setkey(frag.LUT, "Reads")

# Remove names from frag.LUT
frag.LUT[, `:=`(c("Names", "Frequency"), NULL)]

# # Merge temp.table.multi with table.blastn using by = .EACHI to avoid large allocation
# temp.table.multi <- temp.table.multi[frag.LUT, nomatch = 0, allow.cartesian = TRUE]

# Merge temp.table.multi with table.blastn using by = .EACHI to avoid large allocation
multi_traced_reads <- multi_traced_reads[frag.LUT, nomatch = 0, allow.cartesian = TRUE]

# # Remove Type from temp.table.multi
# temp.table.multi[, `:=`(Type, NULL)]

# Remove Type from temp.table.multi
multi_traced_reads[, `:=`(Type, NULL)]

# setkeyv(temp.table.multi, c("BC", "LUTnr"))

setkeyv(multi_traced_reads, c("BC", "LUTnr"))

# # Aggregate column by key
# temp.table.multi[, `:=`(mCount, sum(mCount)), by = key(temp.table.multi)]

# Aggregate column by key
multi_traced_reads[, `:=`(mCount, sum(mCount)), by = key(multi_traced_reads)]

# # Remove duplicates from temp.table.multi
# temp.table.multi <- unique(temp.table.multi, by = c("BC", "LUTnr"))

# Remove duplicates from temp.table.multi
multi_traced_reads <- unique(multi_traced_reads, by = c("BC", "LUTnr"))

# setkeyv(temp.table.multi, "BC")

setkeyv(multi_traced_reads, "BC")

# # Select rows with the highest mCount
# temp.table.multi <- temp.table.multi[temp.table.multi[, .I[mCount == max(mCount)],
#                                                       by = key(temp.table.multi)]$V1]

# Select rows with the highest mCount
multi_traced_reads <- multi_traced_reads[multi_traced_reads[, .I[mCount == max(mCount)],
                                                      by = key(multi_traced_reads)]$V1]

# # Set Mode to "Amb" for mCount = 1
# temp.table.multi[temp.table.multi$mCount == 1, Mode := "Amb"]

# Set Mode to "Amb" for mCount = 1
multi_traced_reads[multi_traced_reads$mCount == 1, Mode := "Amb"]

# print(paste("Number of barcodes with false mCount:", nrow(temp.table.multi[mCount > tCount])))

print(paste("Number of barcodes with false mCount:", nrow(multi_traced_reads[mCount > tCount])))

# # Combine temp.table.multi and temp.table.multi.clean
# temp.table.multi.consensus <- rbind(temp.table.multi, temp.table.multi.clean, fill = TRUE)
# # Combine temp.table.multi and temp.table.multi.clean
# temp.table.multi.consensus <- bind_rows(temp.table.multi, temp.table.multi.clean)

# Combine temp.table.multi and temp.table.multi.clean
temp.table.multi.consensus <- bind_rows(multi_traced_reads, traced_reads)

print(paste("Total number of definitive Barcodes:", length(grep("Def", temp.table.multi.consensus$Mode))))
print(paste("Total number of ambiguous Barcodes:", length(grep("Amb", temp.table.multi.consensus$Mode))))
# print(paste("Total number of single-read Barcodes:", nrow(temp.table.single)))
print(paste("Total number of single-read Barcodes:", nrow(single_traced_reads)))

# # Combine all data tables
# output.Table <- bind_rows(temp.table.multi.consensus, temp.table.single)

# Combine all data tables
output.Table <- bind_rows(temp.table.multi.consensus, single_traced_reads)

# Count the occurrences of each read
frag_counts <- output.Table[, .N, by = Reads]
setnames(frag_counts, "N", "ReadCount")

# Join the read counts to the output.Table
output.Table <- merge(output.Table, frag_counts, by = "Reads", all.x = TRUE)

# Count the occurrences of each barcode for each fragment
barcode_counts <- output.Table[, .N, by = .(Reads, BC)]
setnames(barcode_counts, "N", "BarcodeCount")

# Identify the most frequent barcode for each fragment
most_frequent_barcodes <- barcode_counts[order(-BarcodeCount), .SD[1], by = Reads]

print("Most frequent barcodes for each fragment:")
print(head(most_frequent_barcodes))

# Assuming output.Table is already defined and contains the combined data

# Identify the most frequent read
most_frequent_reads <- frag_counts[order(-ReadCount), .SD[1]]

print("Most frequent reads:")
print(head(most_frequent_reads))

# Save the output table as CSV
write.csv(output.Table, file = "02_analysis/f_multipleContfragmentsComplete.csv", row.names = FALSE)

# Zipping the output table
file_to_zip <- "02_analysis/f_multipleContfragmentsComplete.csv"

# Define the name of the zip file
zip_file <- "02_analysis/f_multipleContfragmentsComplete.zip"

# Zip the file
zip(zip_file, file_to_zip)

# Save the output table as RDA
save(output.Table, file = "02_analysis/f_multipleContfragmentsComplete.rda")

