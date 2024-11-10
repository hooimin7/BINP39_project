#' ---
#' title: "Generate a complete library range object"
#' author: "Tomas Bjorklund"
#' edited by: "Hooi Min Tan Grahn"
rm(list = ls())

library(data.table)
library(ShortRead)
library(digest)
library(scales)
library(Biostrings)
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(devtools))


#'Load sequences
#'===================
# Read LUT.dna CSV file
LUT.dna <- fread("LUTdna.csv")

load("../05_bowtie/alignedLibraries.rda")

load("02_analysis/multipleContfragmentsComplete.rda")

setkey(output.Table, LUTnr)

setkey(LUT.dna, LUTnr)

# Merge the two tables
output.Table <- output.Table[LUT.dna, nomatch = 0]

# Remove unnecessary columns
output.Table[, `:=`(c("Names", "i.NewOrOld"), NULL)]

setnames(output.Table, "Seq", "fragment")

setkey(output.Table, fragment)

# create a new data.table with a idxFrag column which is the index of the fragments
# when there are multiple fragments due to multiple samples taken
range.idx <- data.table(fragment=mcols(allFragments.ranges)$Sequence, 
                        idxFrag=1:length(allFragments.ranges), key="fragment")
# only keep the fragments that are in the range.idx                         
output.Table <- output.Table[range.idx, nomatch=0, allow.cartesian=TRUE]

# create a new data.table with the fragment that are found
foundFragments.ranges <- allFragments.ranges[output.Table$idxFrag]
# remove unnecessary columns
output.Table[,c("Reads","fragment","idxFrag","LUTnr", "Type"):=NULL]
# create a new column RNAcount with the count of the fragment
output.Table[,RNAcount:=tCount]

# add the output.Table to the foundFragments.ranges
mcols(foundFragments.ranges) <- c(mcols(foundFragments.ranges), output.Table)

saveRDS(foundFragments.ranges, file = "02_analysis/completeLibraryRanges_p005.rds")

# Save the merged and annotated sequences to a csv file
write.csv(foundFragments.ranges, file = "02_analysis/completeLibraryRanges_p005.csv", row.names = FALSE)
