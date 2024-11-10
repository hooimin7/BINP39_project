# Clear the workspace
rm(list = ls())

# Increase the maximum allowed size for globals
options(future.globals.maxSize = 8 * 1024^3)  # 8 GiB

# Load necessary libraries
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(future.apply))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(BiocGenerics))
suppressPackageStartupMessages(library(GenomicRanges))

# Set the future plan
plan(multicore, workers = 8)

# Base directory
base_dir <- "/home/hooimin/lu2024-17-19/Data_RAAV/p005_data"

# Define the path to file1
file1 <- file.path(base_dir, "02_rds", "found.p005_72_S59_S8_L001_R1_001.rds")

# Load file1
file1_data <- readRDS(file1)

# Change the default behavior to induce start codons and Methionine
GENETIC_CODE_ALT <- GENETIC_CODE
attr(GENETIC_CODE_ALT, "alt_init_codons") <- c("TAA", "TAG")

# Define the translation function
translate_function <- function(sequence) {
  as.character(Biostrings::translate(DNAString(sequence), genetic.code = GENETIC_CODE_ALT, if.fuzzy.codon = "solve"))
}

# Initialize readCounts
readCounts <- data.table(Group = character(), GroupCount = numeric())

# Function to normalize file1 against a chunk of file2 and file3
normalize_chunk <- function(chunk_index, file1_data, GENETIC_CODE_ALT, readCounts) {
  file2_chunk <- file.path(base_dir, "02_rds", "01_chunkLib", paste0("output_chunks_", chunk_index, ".rds"))
  file3_chunk <- file.path(base_dir, "02_rds", "02_chunkAAV", paste0("found.p005_AAV_01_chunk_", chunk_index, ".rds"))
  
  # Load the chunks
  file2_data <- readRDS(file2_chunk)
  file3_data <- readRDS(file3_chunk)
  
  # Combine the data
  combined_data <- list(file1_data, file2_data, file3_data)
  
  # Normalize the combined data
  # (Insert your normalization logic here)
  
  # Combine the file paths into a character vector
  in.names.all <- c(file1, file2_chunk, file3_chunk)

  # Clean the names
  cleaned_names <- gsub("-", "_", gsub("found\\.|\\.rds", "", basename(in.names.all)))

  # Create a data.table with two columns: Name and GroupName
  cleaned_data <- data.table(Name = cleaned_names, GroupName = "")

  # Define the group name patterns and assign group names
  cleaned_data[grepl("output_chunks_", Name), GroupName := "library"]
  cleaned_data[grepl("_72_", Name), GroupName := "RAAV-72"]
  cleaned_data[grepl("_AAV", Name), GroupName := "AAV"]

  # Assign cleaned_data to load.list 
  load.list <- cleaned_data

  # Assign column names
  setnames(load.list, c("Name", "GroupName"))

  # Select the files to load
  select.Cases <- c(unlist(sapply(load.list$Name, function(x) grep(x, in.names.all), simplify = TRUE)))
  select.Cases <- unique(select.Cases)

  # Load the files
  if (length(select.Cases) > 0) {
    in.names.all <- in.names.all[select.Cases]
    # Ensure select.Cases has names
    names(select.Cases) <- load.list$Name[select.Cases]
  } else {
    in.names.all <- NULL
  }

  # Load the files
  in.names.all <- in.names.all[select.Cases]

  # Ensure select.Cases has names
  names(select.Cases) <- load.list$Name[select.Cases]

  # Create the grouping data frame
  grouping <- data.frame(
    Sample = gsub("-", "_", gsub("found.|(02_rds/)|(.rds)", "", basename(in.names.all))),
    Group = load.list$GroupName[match(names(select.Cases), load.list$Name)],
    stringsAsFactors = FALSE
  )

  # Function to load the RDS files and handle different formats
  loadRDS <- function(in.name) {
    this.sample <- readRDS(in.name)
    
    this.name <- gsub("-", "_", gsub("found.|(02_rds/)|(.rds)", "", basename(in.name)))
    this.group <- grouping[match(this.name, grouping$Sample), "Group"]
    
    if (inherits(this.sample, "GAlignments")) {
      # Handle GAlignments object
      mcols(this.sample) <- cbind(mcols(this.sample),
                                  data.frame(Sample = this.name, Group = this.group,
                                             stringsAsFactors = FALSE))
    } else if (inherits(this.sample, c("data.table", "data.frame"))) {
      # Handle data.table or data.frame object
      this.sample[, Sample := this.name]
      this.sample[, Group := this.group]
    } else {
      warning(paste("Skipping file:", in.name, "- Unsupported object type"))
      return(NULL)
    }
    
    return(this.sample)
  }

  # Load all RDS files specified in 'in.names.all' into a list of R objects, filtering out NULL results
  all.samples <- lapply(in.names.all, loadRDS)
  all.samples <- Filter(Negate(is.null), all.samples)

  # Combine the list of R objects into a single GAlignmentsList object. Using unlist() to flatten the list of objects before combining
  all.samples <- GAlignmentsList(unlist(all.samples))

  # Rename the samples into unique and R valid names
  names(all.samples) <- make.names(names(all.samples), unique = TRUE)

  # Create a data.table with the seqnames and seq length and set the key to seqnames  
  length.Table <- data.table(seqnames = names(seqlengths(all.samples)), seqlength = seqlengths(all.samples), key = "seqnames")

  # Create a data.table from all samples and set the key to seqnames    
  all.samples <- data.table(as.data.frame(all.samples), key = "seqnames")

  # Remove unnecessary columns
  all.samples[, `:=`(c("strand", "qwidth", "cigar", "njunc", "end", "Frequency", "i.Name"), NULL)]

  # Merge all samples with the length table to get the sequence length of every sequence
  all.samples <- all.samples[length.Table] # A data.table merge to match seqlengths to their respective seqnames

  # Create new columns and assign the values based on the seqnames that was split into its parts
  all.samples[, `:=`(c("Category", "Protein", "Origin", "Extra", "Number", "GeneName"), tstrsplit(seqnames, ",", fixed = TRUE))]

  # Remove unnecessary columns
  all.samples[, `:=`(c("seqnames", "Protein", "Origin", "Extra", "Number"), NULL)]

  # Change the GeneName column to hold the actual gene name
  all.samples[, `:=`(GeneName, gsub("#\\d+", "", GeneName))]

  # Save the normalized data
  saveRDS(all.samples, file = paste0("03_normalize/S59_S8/allSamples_p005_S59_chunk", chunk_index, ".RDS"))

  # Load the final combined normalized data
  final_normalized_data <- all.samples

  print("initial")
  print(head(final_normalized_data))

  # Set the key to the group
  setkey(final_normalized_data, Group)

  # Get rid of single read samples
  final_normalized_data <- final_normalized_data[RNAcount > 1, ]

  # Store the read counts grouped by the group
  chunk_readCounts <- final_normalized_data[, list(GroupCount = sum(RNAcount)), by = "Group"]

  # Update readCounts with the current chunk
  readCounts <- rbind(readCounts, chunk_readCounts, fill = TRUE)

  # Normalize the group reads by the max Group reads set values between 0 and 1 where 1 is the max value
  readCounts[, GroupCount := GroupCount / max(GroupCount)]

  # Save intermediate readCounts
  saveRDS(readCounts, file = paste0("03_normalize/S59_S8/readCounts_p005_S59_chunk", chunk_index, ".RDS"))

  # Normalize the read counts by the Group Count to correct for variable read depth
  final_normalized_data <- final_normalized_data[readCounts]
  final_normalized_data[, RNAcount := RNAcount / GroupCount]

  # Extract only the definitive Barcodes which have at least 2 reads per pair of barcode & LUTnr
  final_normalized_data <- final_normalized_data[Mode == "Def"]

  # Save intermediate normalized data
  saveRDS(final_normalized_data, file = paste0("03_normalize/S59_S8/normalized_p005_S59_data", chunk_index, ".RDS"))

  # Set the key to the Group
  setkey(final_normalized_data, Group)

  # Get all the AAV samples by removing the DNA samples
  total.AAV.samples <- final_normalized_data[Group != "library"]

  # Get the transported AAV samples
  transported.AAV.samples <- total.AAV.samples[grepl("AAV", total.AAV.samples$Group)]

  # Get the transported 72 tissue AAV samples
  transported.AAV.samples.72 <- total.AAV.samples[grepl("RAAV-72",total.AAV.samples$Group)]
  # Set the group name to mRNA_72
  transported.AAV.samples.72[,Group := "mRNA_72"]

  # Set the group name to mRNA_All
  total.AAV.samples[, Group := "mRNA_All"]

  print("total.AAV.samples")
  print(head(total.AAV.samples))

  # Recombine all samples with their new group names
  final_normalized_data <- rbind(final_normalized_data, total.AAV.samples, transported.AAV.samples.72)

  print("rbind")
  print(head(final_normalized_data))

  # Remove the data.tables that were created before
  rm(total.AAV.samples, transported.AAV.samples.72)

  # Set the key of all samples to a vector of the columns
  setkeyv(final_normalized_data, c("Group", "Category", "GeneName", "start", "width", "Sequence", "seqlength"))

  # Now combine all information about identical fragments that are found with different barcodes
  final_normalized_data <- final_normalized_data[, j = list(mCount = sum(mCount),
                                                            tCount = sum(tCount),
                                                            BC = paste(unique(BC), collapse = ","),
                                                            LUTnrs = paste(unique(LUTnr), collapse = ","),
                                                            RNAcount = sum(RNAcount),
                                                            NormCount = log2(sum(RNAcount) + 1) * .N),
                                                 by = c("Group", "Category", "GeneName", "start", "width", "Sequence", "seqlength")]

  print("identical")
  print(head(final_normalized_data))

  # Adjust the start, width, and seqlength from DNA based to AA based
  final_normalized_data[, start := floor((start + 2) / 3)]
  final_normalized_data[, width := ceiling((width) / 3)]
  final_normalized_data[, seqlength := ceiling(seqlength / 3)]

  # Calculate the absolute AA position of the middle of the fragment 
  final_normalized_data[, AA := floor(start + (width / 2))]

  # Calculate the relative AA position of the middle of the fragment in the whole sequence
  final_normalized_data[, AAproc := AA / seqlength * 100]

  print("AAproc")
  print(head(final_normalized_data))

  # Apply translation function in parallel
  final_normalized_data[, Peptide := future_lapply(Sequence, translate_function)]

  # Change the Peptide column to character
  final_normalized_data[, Peptide := as.character(Peptide)]

  print("final")
  print(head(final_normalized_data))

  # Print the number of samples in the final normalized data table
  print(nrow(final_normalized_data))

  # Save the final normalized data
  saveRDS(final_normalized_data, file = paste0("03_normalize/S59_S8/allSamplesDataTable_p005_72_S59_chunk", chunk_index, ".RDS"))
}

# Process each chunk
num_chunks <- 39  # Total number of chunks

for (i in 1:num_chunks) {
  normalize_chunk(i, file1_data, GENETIC_CODE_ALT, readCounts)
}

for (i in 1:num_chunks) {
  file.remove(paste0("03_normalize/S59_S8/allSamples_p005_S59_chunk", i, ".RDS"))
  file.remove(paste0("03_normalize/S59_S8/readCounts_p005_S59_chunk", i, ".RDS"))
  file.remove(paste0("03_normalize/S59_S8/normalized_p005_S59_data", i, ".RDS"))
}