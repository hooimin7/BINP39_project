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

# Start time
strt1 <- Sys.time()

# Define the base directory
base_dir <- "/home/hooimin/lu2024-17-19/Data_RAAV/p005_data"

LUT.dna <- fread("../../RAAV-60/p005/04_blast/LUTdna.csv")

# Define the path to the 03_normalize subdirectory
normalize_dir <- file.path(base_dir, "03_normalize")
  
# List all subdirectories within the 03_normalize directory
subdirectories <- list.dirs(path = normalize_dir, full.names = TRUE, recursive = TRUE)

# Filter out the 03_normalize directory itself
subdirectories <- subdirectories[subdirectories != normalize_dir]

# Print the result
print(subdirectories)
  
# Loop through each subdirectory
for (subdir in subdirectories) {
  # Extract the base name of the subdirectory
  base_name <- basename(subdir)

  # List RDS files matching the pattern in the current subdirectory
  rds.files <- list.files(
    path = subdir,
    pattern = "allSamplesDataTable_p005_(26_S3|27_S4|28_S5|30_S6|31_S7|32_S8|34_S9|35_S10|36_S11|37_S12|38_S13|39_S14|42_S15|43_S16|44_S17|46_S18|49_S19|50_S20|51_S21|52_S22)\\.RDS$",
    full.names = TRUE
  )

  # Check if there is exactly one RDS file
  if (length(rds.files) == 1) {
    # Read the single RDS file
    rds_files <- readRDS(rds.files[1])

    setkey(rds_files, Group)

    seq.arry <- LUT.dna$LUTnr

    print("seq.arry")
    print(head(seq.arry))

    # Extract unique LUTnr values from the library group
    seq.lib <- unique(rds_files[J("library")]$LUTnr)

    # Split concatenated LUTnr values in seq.lib into individual values
    seq.lib_split <- unlist(strsplit(seq.lib, ","))

    # Save seq.lib_split
    saveRDS(seq.lib_split, file = file.path(subdir, paste0("seq_lib_split_", base_name, ".RDS")))

    # Extract unique LUTnr values from the AAV group
    seq.AAV <- unique(rds_files[J("AAV")]$LUTnr)

    # Save seq.AAV
    saveRDS(seq.AAV, file = file.path(subdir, paste0("seq_AAV_", base_name, ".RDS")))

    # Extract unique LUTnr values from the mRNA_53 group
    seq.53 <- unique(rds_files[grep("mRNA_53", Group)]$LUTnr)

    # Save seq.53
    saveRDS(seq.53, file = file.path(subdir, paste0("seq_53_", base_name, ".RDS")))

  }
}