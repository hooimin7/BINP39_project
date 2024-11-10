#' ---
#' title: "Pairwise sample analysis output"
#' author: "Tomas Bjorklund"
#' edited by: "Hooi Min"
#' output: html_document
#' ---

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(devtools))

# Define the base directory
base_dir <- "/home/hooimin/lu2024-17-19/Data_RAAV/p007_data"

# Define the path to the 06_motif subdirectory
motif_dir <- file.path(base_dir, "06_motif")
  
# List all subdirectories within the 06_motif directory
subdirectories <- list.dirs(path = motif_dir, full.names = TRUE, recursive = TRUE)

# Filter out the 06_motif directory itself
subdirectories <- subdirectories[subdirectories != motif_dir]

# Initialize an empty list to store data from all subdirectories
data_list <- list()

# Loop through each subdirectory
for (subdir in subdirectories) {
  # Extract the base name of the subdirectory
  base_name <- basename(subdir)

  # Define the path to the Hammock subdirectory
  hammock_dir <- file.path(subdir, "Hammock")

  # List TSV files matching the pattern in the Hammock subdirectory
  tsv_files <- list.files(
    path = hammock_dir,
    pattern = "final_clusters.tsv",
    full.names = TRUE
  )

  # Check if there is exactly one TSV file
  if (length(tsv_files) == 1) {
    motif_clusters_all_path <- tsv_files[1]

    # Read the TSV file into a data.table
    motif.all <- data.table(read.table(motif_clusters_all_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE))

    # Set the key to 'sum'
    setkey(motif.all, sum)

    # Select the top 25 samples
    motif.all.top <- motif.all[, head(.SD, 25), by = sum]

    # Add the top 25 samples to the data list
    data_list[[length(data_list) + 1]] <- motif.all.top

    # Save the top 25 samples to a TSV file
    write.table(motif.all.top, file = file.path(subdir, paste0("25_MostMotif_", base_name, ".tsv")), sep = "\t", row.names = FALSE, col.names = TRUE)
  }
}

# Concatenate all data frames from different subdirectories into one
concatenated_data <- rbindlist(data_list, fill = TRUE)

# Sort by 'sum' in descending order
concatenated_data <- concatenated_data[order(-sum)]

# Save the top 25 samples to a TSV file
write.table(concatenated_data, file = file.path(motif_dir, "top25_p007_53.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE)