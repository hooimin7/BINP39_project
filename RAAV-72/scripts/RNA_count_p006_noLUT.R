# Clear the workspace
rm(list = ls())

# Set the working directory to the base directory
setwd("/home/hooimin/lu2024-17-19/RAAV-72")

# Load necessary libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(future.apply))

# Increase the maximum allowed size for globals
options(future.globals.maxSize = 20 * 1024^3)  # 20 GiB

# Base directory
base_dir <- "/home/hooimin/lu2024-17-19/RAAV-72"

# List of directories ending with _R1_001
directories <- list.dirs(base_dir, recursive = TRUE, full.names = TRUE)
directories <- directories[grepl("_R1_001$", directories)]

# Function to process each directory
process_directory <- function(dir) {
  # Extract the base name of the directory
  dir_base <- basename(dir)
  out_dir <- file.path(dir, "06_traceBack", "02_analysis")

  # Ensure the output directory exists
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # File paths for barcodes and starcodes
  barcodes_file <- file.path(dir, "03_pairfq", "barcodes_6_BarLib_paired_reads.fastq.gz")
  starcode_file <- file.path(dir, "05_starcode", "01_data", "barcodes_6_reduced.txt")

  # Check if files exist
  if (file.exists(barcodes_file) && file.exists(starcode_file)) {
    # Load necessary data
    load("../RAAV-60/p006/04_blast/02_analysis/multipleContfragmentsCompleteNoLUT.rda")

    # Read FASTQ files
    reads.BC <- readFastq(barcodes_file)

    barcodeTable <- data.table(ID = as.character(ShortRead::id(reads.BC)), BC = as.character(sread(reads.BC)), key = "BC")

    # Reduce barcode
    table.BC.sc <- data.table(read.table(starcode_file, header = FALSE, row.names = 1, skip = 0, sep = "\t", stringsAsFactors = FALSE, fill = FALSE), keep.rownames = TRUE, key = "rn")
    table.BC.sc[, `:=`(V2, NULL)]
    table.BC.sc <- table.BC.sc[, strsplit(as.character(V3), ",", fixed = TRUE), by = rn]

    unique_reads_BC <- length(unique(sread(reads.BC)))
    unique_table_BC_sc <- length(unique(table.BC.sc$V1) %in% unique(sread(reads.BC)))
    SC.droppedBC <- unique_reads_BC - unique_table_BC_sc

    print(paste("SC.droppedBC:", SC.droppedBC))

    setnames(table.BC.sc, c("V1", "rn"), c("BC", "scBC"))

    # Replacing barcodes with Starcode reduced versions
    setkey(table.BC.sc, BC)
    setkey(barcodeTable, BC)

    # Merging data tables
    barcodeTable <- barcodeTable[table.BC.sc, nomatch = 0]

    # Renaming columns
    setnames(barcodeTable, c("BC", "scBC"), c("oldBC", "BC"))
    setkey(barcodeTable, BC)

    # Calculating unique barcodes
    allBCs <- length(unique(barcodeTable$oldBC))
    print(paste("allBCs, Original unique barcodes:", dir_base, allBCs))

    scBC <- length(unique(barcodeTable$BC))
    print(paste("scBC, SC reduced unique barcodes:", dir_base, scBC))

    invisible(barcodeTable[, `:=`(oldBC, NULL)])

    setkey(output.Table, "BC")

    # Creating and formatting tables
    BCcount <- data.table(as.data.frame(rev(sort(table(barcodeTable$BC))), row.names = "Var1"), keep.rownames = TRUE)
    setnames(BCcount, colnames(BCcount), c("BC", "RNAcount"))
    setkey(BCcount, "BC")

    foundFrags <- output.Table[BCcount, nomatch = 0]

    print("foundFrags")
    print(head(foundFrags))

    print("outtable")
    print(head(output.Table))

    # Check the column names in foundFrags and output.Table
    print("Column names in foundFrags:")
    print(colnames(foundFrags))

    print("Column names in output.Table:")
    print(colnames(output.Table))

    # Function to index fragments that are found multiple times since they are the same in multiple sequences from the known Retrograde_transport seq
    matchRange <- function(idxFrag) {
      matchRanges <- which(output.Table$Reads == foundFrags$Reads[idxFrag])
      if (length(matchRanges) == 0) {
        return(matrix(ncol = 2, nrow = 0))  # Return an empty matrix if no matches
      }
      result <- cbind(matchRanges, idxFrag)
      result <- as.matrix(result)  # Ensure the result is a matrix 
      return(result)
    }

    # Setup parallel backend using future
    plan(multisession, workers = detectCores() / 2)

    # Create a list of the matches
    match.ranges.list <- future_lapply(1:nrow(foundFrags), function(i) {
        matchRange(i)
    })

    # Ensure match.ranges.list is a list
    if (!is.list(match.ranges.list)) {
      match.ranges.list <- list(match.ranges.list)
    }

    # Create a matrix of the matches
    match.ranges <- do.call(rbind, match.ranges.list)

    # Check if match.ranges is not empty and has the correct structure
    if (nrow(match.ranges) == 0 || ncol(match.ranges) != 2) {
      print("No valid matches found for this chunk.")
      next
    }

    # Create a list of the found fragments
    foundFragments.ranges <- output.Table[match.ranges[, 1]]

    # Define the output name
    name.out <- paste0('p006_72_', dir_base)

    # If there are more than one match, then save the found fragments
    if (ncol(match.ranges) >= 2) {
      foundFrags <- foundFrags[match.ranges[, "idxFrag"], ]

      # Remove unnecessary columns
      foundFrags[, c("Reads", "fragment", "LUTnr") := NULL]

      # Add found fragments to the foundFragments.ranges
      foundFragments.ranges <- cbind(foundFragments.ranges, foundFrags)
      # Sort them by the RNA count in descending order
      foundFragments.ranges <- foundFragments.ranges[order(-foundFragments.ranges$RNAcount)]

      print("sort RNA")
      print(head(foundFragments.ranges))

      # Save the found fragments for the sample in the output folder
      saveRDS(foundFragments.ranges, file = file.path(out_dir, paste0("found.", name.out, ".rds")), compress = TRUE)

      # Save the foundFrags as CSV
      write.csv(foundFragments.ranges, file = file.path(out_dir, paste0("found.", name.out, ".csv")), row.names = FALSE)
    } else {
      message("Files not found in directory: ", dir)
    }
  }
}

# Run RNA count on each directory
lapply(directories, process_directory)