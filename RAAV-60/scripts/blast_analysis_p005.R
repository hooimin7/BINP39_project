#' ---
#' title: "Blast analysis"
#' author: "Hooi Min Tan Grahn"
rm(list = ls())

library(data.table)
library(utils)
library(ShortRead)
library(digest)
library(scales)
library(Biostrings)

#'Load sequences
#'===================
# Read LUT.dna CSV file
LUT.dna <- fread("LUTdna.csv")

# Create ShortRead object
LUT.seq <- ShortRead(DNAStringSet(LUT.dna$Seq), BStringSet(1:length(LUT.dna$LUTnr)))

## Save LUT.seq to a file 
# writeXStringSet(LUT.seq@sread, filepath = "02_analysis/LUT_seq.fasta")

# File paths for fragments and barcodes
fragments.file <- "../06_library/fragments_5_FragBar_paired_again.fastq.gz"
barcodes.file <- "../06_library/barcodes_5_FragBar_paired_again.fastq.gz"

# Read FASTQ files
reads.trim <- readFastq(fragments.file)
reads.BC <- readFastq(barcodes.file)

unique.reads <- unique(sread(reads.trim))
unique.reads <- ShortRead(DNAStringSet(unique.reads), BStringSet(1:length(unique.reads)))

## Save unique reads to a file
# writeXStringSet(unique.reads@sread, filepath = "02_analysis/unique_reads.fasta")

print("unique.reads")
print(head(unique.reads))

# Read blast output into data.table
table.blastn <- data.table(scan(file = "01_blastout/blast_output_5.csv.gz", what = "character",
sep = ";"), keep.rownames = FALSE, key = "V1")

print("Initial table_blastn:")
print(head(table.blastn))

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

print("After splitting columns:")
print(head(table.blastn))

# Create lookup tables for Reads and Sequence
unique_reads <- unique(table.blastn$Reads)
unique_sequences <- unique(table.blastn$Sequence)

reads_lookup <- setNames(seq_along(unique_reads), unique_reads)
sequence_lookup <- setNames(seq_along(unique_sequences), unique_sequences)

# Map original values to integer indices
table.blastn[, Reads := reads_lookup[Reads]]
table.blastn[, Sequence := sequence_lookup[Sequence]]

# Update Reads and Sequence columns with the corresponding sequences
table.blastn[, Reads := as.character(sread(unique.reads[as.integer(Reads)]))]
table.blastn[, Sequence := as.character(sread(LUT.seq[as.integer(Sequence)]))]

print("After updating Reads and Sequence:")
print(head(table.blastn))


# Set keys for merging
setkey(table.blastn, Sequence)
setkey(LUT.dna, Seq)

# Merge table.blastn with LUT.dna using .EACHI
table.blastn <- table.blastn[LUT.dna, nomatch = 0, allow.cartesian = TRUE]

print("After merging:")
print(head(table.blastn))

# Remove unnecessary columns
table.blastn[, `:=`(c("V1", "identity", "alignmentLength", "gapOpens", "q_start",
"q_end", "s_start", "s_end", "evalue", "Sequence", "Names"), NULL)]

# Convert bitScore and mismatches to numeric
table.blastn[, `:=`(bitScore, as.numeric(bitScore))]
table.blastn[, `:=`(mismatches, as.numeric(mismatches))]

print("After converting bitScore and mismatches to numeric:")
print(head(table.blastn))

# Set keys and order for unique alignment
setkeyv(table.blastn, c("Reads", "LUTnr"))
setorder(table.blastn, Reads, LUTnr, -bitScore) #This makes sure that a fragment is only aligned once to the table.blastn <- unique(table.blastn, by = c("Reads", "LUTnr"))

# Ensure unique alignment
table.blastn <- unique(table.blastn, by = c("Reads", "LUTnr"))

print("After ensuring unique alignment:")
print(table.blastn)

# Get top hit alignment
table.blastn.topHit <- table.blastn[table.blastn[, .I[which.max(bitScore)],
by = "Reads"]$V1]

print("tophit")
print(table.blastn.topHit)

# Create full table with reads and barcodes
full.table <- data.table(Reads = as.character(sread(reads.trim)), BC = as.character(sread(reads.BC)),
key = "Reads")


print("full_table")
print(full.table)

# Get total number of reads
all.reads <- nrow(full.table)

print("all.reads")
print(all.reads)

# Merge full.table with top hit alignment
full.table <- full.table[table.blastn.topHit, nomatch = 0] # Merge reads with the top hit alignment

print("full_table_tophit")
print(full.table)

# Print alignment percentage
print(paste("Alignment percentage:", percent(nrow(full.table)/all.reads)))

# Reduce barcode
out.name.BC.star <- "../07_starcode/01_data/barcodes_5_reduced.txt"

table.BC.sc <- data.table(read.table(out.name.BC.star, header = FALSE, row.names = 1,
skip = 0, sep = "\t", stringsAsFactors = FALSE, fill = FALSE), keep.rownames = TRUE,
key = "rn") #, nrows = 1000

print("initial table BC ")
print(table.BC.sc)

table.BC.sc[, `:=`(V2, NULL)]

table.BC.sc <- table.BC.sc[, strsplit(as.character(V3), ",", fixed = TRUE),
by = rn]

print("V3, rn table BC ")
print(table.BC.sc)

# SC.droppedBC <- length(unique(sread(reads.BC))) - length(unique(table.BC.sc$V1) %in% unique(sread(reads.BC)))

unique_reads_BC <- length(unique(sread(reads.BC)))
unique_table_BC_sc <- length(unique(table.BC.sc$V1) %in% unique(sread(reads.BC)))
SC.droppedBC <- unique_reads_BC - unique_table_BC_sc

print(unique_reads_BC)
print(unique_table_BC_sc)
print(paste("Dropped BCs in Starcode:", SC.droppedBC))

setnames(table.BC.sc, c("V1", "rn"), c("BC", "scBC"))


print("Column names in table.BC.sc:")
print(colnames(table.BC.sc))

# Replacing barcodes with Starcode reduced versions
# Setting keys for data tables
setkey(full.table, BC)
setkey(table.BC.sc, BC)

# Merging data tables
full.table <- full.table[table.BC.sc, nomatch = 0]
# full.table <- merge(full.table,table.BC.sc, by='BC', all = FALSE, all.x =
# FALSE) rm(table.BC.sc)
print("full.table merged with BC.sc")
print(full.table)

# Renaming columns
setnames(full.table, c("BC", "scBC"), c("oldBC", "BC"))
setkey(full.table, BC)

# Calculating unique barcodes
RetainedBC <- length(unique(full.table$oldBC))
scBC <- length(unique(full.table$BC))
print(paste("Original unique barcodes:", RetainedBC))
print(paste("SC reduced unique barcodes:", scBC))

# Creating and formating tables
table.frag <- data.table(as.data.frame((rev(sort(table(full.table$oldBC))))[1:10],
row.names = "Var1"), keep.rownames = TRUE)
print("table.frag")
print(table.frag)
# In R versions below 3.3 remove, row.names = 'Var1' to make this compatible

setnames(table.frag, colnames(table.frag), c("Original BC", "Count"))

table.frag <- data.table(as.data.frame((rev(sort(table(full.table$BC))))[1:10],
row.names = "Var1"), keep.rownames = TRUE)
# In R versions below 3.3 remove, row.names = 'Var1' to make this compatible
setnames(table.frag, colnames(table.frag), c("SC reduced BC", "Count"))

# Removing oldBC column
invisible(full.table[, `:=`(oldBC, NULL)])

print("table.frag_after sc")
print(table.frag)
print(full.table)


# Ordering the table
full.table <- full.table[order(full.table$BC), ]

print("full.table ordering")
print(full.table)

# Converting mismatches to numeric
full.table[, `:=`(mismatches, as.numeric(mismatches))]

# Splitting reads into single-read and multi-read barcodes
temp.table.single <- full.table[full.table[, .I[.N == 1], by = "BC"]$V1] # BC appears once 
temp.table.multi <- full.table[full.table[, .I[.N > 1], by = "BC"]$V1] # BC appears more than once

print("spliting")
print(temp.table.single)
print(temp.table.multi)

# Adding columns to single read table
temp.table.single[, `:=`(c("mCount", "tCount"), 1)] # Add mCount and tCount col with a vaule of 1
temp.table.single$Mode <- "Amb" # Add a Mode column with value (Ambiguos)

print("after adding col to single")
print(temp.table.single)

# Setting key for multireads table
setkeyv(temp.table.multi, c("BC", "LUTnr"))

# Aggregating multi-read table
temp.table.multi[, `:=`(c("bitScore", "mismatches", "tCount"), list(mean(bitScore),
median(mismatches), .N)), by = key(temp.table.multi)] # Set tCount to the number of rows in each group
temp.table.multi$Mode <- "Def" # Add definitive column
temp.table.multi <- unique(temp.table.multi) # Remove duplicate rows

print("after aggregating multi-read table")
print(temp.table.multi)


print("Utilized reads.......")
print(nrow(full.table))
print("Whereof single reads.......")
print(nrow(temp.table.single))

# Splitting multi-read barcodes into clean and chimeric

setkeyv(temp.table.multi, "BC")

# Filter for barcodes
temp.table.multi.clean <- temp.table.multi[temp.table.multi[, .I[.N == 1], by = "BC"]$V1] # BC appears only once

print("table multi clean")
print(temp.table.multi.clean)

# Filter for non-unique barcodes
temp.table.multi <- temp.table.multi[temp.table.multi[, .I[.N > 1], by = "BC"]$V1] # BC appears more than once

# Update mCount column
temp.table.multi.clean[, `:=`(mCount, tCount)]

print("Clean multi-read barcodes.......")
print(nrow(temp.table.multi.clean))
print("Chimeric multi-read barcodes.......")
print(length(unique(temp.table.multi$BC)))

# Calculate consensus alignment of chimeric barcodes
setkey(temp.table.multi, "BC")

# Update mCount column
temp.table.multi[, mCount := tCount]

print("update mCount")
print(temp.table.multi)

# Aggregate tCount by BC
temp.table.multi[, tCount := sum(tCount), by = "BC"]

print("Aggregate tCount")
print(temp.table.multi)

setkey(temp.table.multi, "Reads")

# Remove columns
temp.table.multi[, `:=`(c("LUTnr", "bitScore", "mismatches", "Type"), NULL)]

print("Remove columns temp.table")
print(temp.table.multi)

# Separate unique and non-unique rows in temp.table.multi
unique_temp_table_multi <- temp.table.multi[!duplicated(temp.table.multi)]
non_unique_temp_table_multi <- temp.table.multi[duplicated(temp.table.multi)]

print("Unique rows in temp.table.multi")
print(unique_temp_table_multi)

print("Non-unique rows in temp.table.multi")
print(non_unique_temp_table_multi)

# Remove columns
table.blastn[, c("Type") := NULL]

print("Remove columns table.blastn")
print(table.blastn)

# Separate unique and non-unique rows in table.blastn
unique_table_blastn <- table.blastn[!duplicated(table.blastn)]
non_unique_table_blastn <- table.blastn[duplicated(table.blastn)]

print("Unique rows in table.blastn")
print(unique_table_blastn)

print("Non-unique rows in table.blastn")
print(non_unique_table_blastn)

setkey(table.blastn, "Reads")



# Merge temp.table.multi with table.blastn using by = .EACHI to avoid large allocation
temp.table.multi <- temp.table.multi[table.blastn, nomatch = 0, allow.cartesian = TRUE, by = .EACHI]

setkeyv(temp.table.multi, c("BC", "LUTnr"))

# Aggregate column by key
temp.table.multi[, `:=`(c("bitScore", "mismatches", "mCount"), list(max(bitScore, na.rm = TRUE),
median(mismatches, na.rm = TRUE), sum(mCount))), by = key(temp.table.multi)]

print("Aggregate col")
print(temp.table.multi)

# Remove duplicates from temp.table.multi
temp.table.multi <- unique(temp.table.multi, by = c("BC", "LUTnr"))

print("Duplicate removal")
print(temp.table.multi)

setkeyv(temp.table.multi, "BC")

# Select rows with the highest mCount
temp.table.multi <- temp.table.multi[temp.table.multi[, .I[mCount == max(mCount)],
by = key(temp.table.multi)]$V1]

print("Highest mCount")
print(temp.table.multi)

# Select only rows with the highest bitScore
temp.table.multi <- temp.table.multi[temp.table.multi[, .I[which.max(bitScore)],
by = key(temp.table.multi)]$V1]

print("Highest bitScore")
print(temp.table.multi)

# Set Mode to "Amb" for mCount = 1
temp.table.multi[temp.table.multi$mCount == 1, Mode := "Amb"]

print(paste("Number of barcodes with false mCount:", nrow(temp.table.multi[mCount > tCount])))

# Combine temp.table.multi and temp.table.multi.clean
temp.table.multi.consensus <- rbind(temp.table.multi, temp.table.multi.clean, fill = TRUE)

print("Table consensus")
print(temp.table.multi.consensus)
print(paste("Total number of definitive Barcodes:", length(grep("Def", temp.table.multi.consensus$Mode))))
print(paste("Total number of ambiguous Barcodes:", length(grep("Amb", temp.table.multi.consensus$Mode))))
print(paste("Total number of single-read Barcodes:", nrow(temp.table.single)))

# Combine all data tables
output.Table <- rbind(temp.table.multi.consensus, temp.table.single, fill = TRUE)

# Save the output table as CSV
write.csv(output.Table, file = "02_analysis/multipleContfragmentsComplete.csv", row.names = FALSE)

# Save the output table as RDA
save(output.Table, file = "02_analysis/multipleContfragmentsComplete.rda")