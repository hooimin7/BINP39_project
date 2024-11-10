rm(list=ls()) # Clear workspace
library(data.table)
library(ShortRead)
library(Biostrings)
library(parallel)
library(doParallel)
library(foreach)
library(iterators)

# Load sequences
# ===================

# Read LUT.dna CSV file
LUT.dna <- fread("../../../RAAV-60/p006/04_blast/LUTdna.csv")

load("../../../RAAV-60/p006/04_blast/02_analysis/multipleContfragmentsComplete.rda")

load("../../../RAAV-60/p005/05_bowtie/alignedLibraries.rda")

# File paths for barcodes
barcodes.file <- "../03_pairfq/barcodes_7_BarLib_paired_reads.fastq.gz"

# Read FASTQ files
reads.BC <- readFastq(barcodes.file)

barcodeTable <- data.table(ID = as.character(ShortRead::id(reads.BC)), BC = as.character(sread(reads.BC)),
key = "BC")

# Reduce barcode
out.name.BC.star <- "../05_starcode/01_data/barcodes_7_reduced.txt"

table.BC.sc <- data.table(read.table(out.name.BC.star, header = FALSE, row.names = 1,
                                     skip = 0, sep = "\t", stringsAsFactors = FALSE, fill = FALSE), keep.rownames = TRUE,
                          key = "rn") #, nrows = 1000

table.BC.sc[, `:=`(V2, NULL)]

table.BC.sc <- table.BC.sc[, strsplit(as.character(V3), ",", fixed = TRUE),
                           by = rn]

unique_reads_BC <- length(unique(sread(reads.BC)))
unique_table_BC_sc <- length(unique(table.BC.sc$V1) %in% unique(sread(reads.BC)))
SC.droppedBC <- unique_reads_BC - unique_table_BC_sc

print(paste("SC.droppedBC:", SC.droppedBC))

setnames(table.BC.sc, c("V1", "rn"), c("BC", "scBC"))

# Replacing barcodes with Starcode reduced versions
# Setting keys for data tables
setkey(table.BC.sc, BC)

setkey(barcodeTable, BC)

# Merging data tables
barcodeTable <- barcodeTable[table.BC.sc, nomatch = 0]

# Renaming columns
setnames(barcodeTable, c("BC", "scBC"), c("oldBC", "BC"))
setkey(barcodeTable, BC)

# Calculating unique barcodes
allBCs <- length(unique(barcodeTable$oldBC))

print(paste("allBCs, Original unique barcodes:", allBCs))

scBC <- length(unique(barcodeTable$BC))

print(paste("scBC, SC reduced unique barcodes:", scBC))

invisible(barcodeTable[, `:=`(oldBC, NULL)])

setkey(output.Table, "BC")

# Creating and formating tables
BCcount <- data.table(as.data.frame(rev(sort(table(barcodeTable$BC))), row.names = "Var1"),
keep.rownames = TRUE)
# In R versions below 3.3 remove, row.names = 'Var1' to make this compatible

setnames(BCcount, colnames(BCcount), c("BC", "RNAcount"))

setkey(BCcount, "BC")

foundFrags <- output.Table[BCcount, nomatch = 0]

setkey(foundFrags, "LUTnr")

setkey(LUT.dna, "LUTnr")

foundFrags <- foundFrags[LUT.dna, nomatch = 0]

setnames(foundFrags, "Seq", "fragment")

foundFrags[, `:=`(c("Names", "i.NewOrOld", "Type"), NULL)]

# Define the output name
name.out <- "p007_AAV_03_frag"

# Save the found fragments for the sample in the output folder
saveRDS(foundFrags, file = paste("02_analysis/", "found.", name.out, ".rds", sep = ""), 
          compress = TRUE)

# Save the foundFrags as CSV
write.csv(foundFragments.ranges, file = paste0("02_analysis/", "found.", name.out, ".csv"), row.names = FALSE)