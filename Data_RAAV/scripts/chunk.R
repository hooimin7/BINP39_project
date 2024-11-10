#' ---
#' title: "Normalize counts"
#' author: "Hooi Min Tan Grahn"
#' ---

# Clear the workspace
rm(list = ls())

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicAlignments))

# Load the large file2
complibFile <- readRDS("02_rds/01_chunkLib/completeLibraryRanges_p007.rds")

# Ensure foundFragments.ranges is defined
if (!exists("foundFragments.ranges")) {
  # Assuming foundFragments.ranges is part of complibFile
  foundFragments.ranges <- complibFile
}

print(head(foundFragments.ranges))

# Split foundFragments.ranges into smaller chunks
chunk_size <- ceiling(length(foundFragments.ranges) / 29)  # Adjust the chunk size to be smaller
output_chunks <- split(foundFragments.ranges, ceiling(seq_along(foundFragments.ranges) / chunk_size))

# Save chunks to disk
for (i in seq_along(output_chunks)) {
  saveRDS(output_chunks[[i]], paste0("02_rds/01_chunkLib/output_chunks_", i, ".rds"))
}

print("Chunks saved successfully.")