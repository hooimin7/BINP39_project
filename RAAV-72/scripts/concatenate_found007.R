# Clear the workspace
rm(list = ls())

# Set the working directory to the base directory
setwd("/home/hooimin/lu2024-17-19/RAAV-72")

# Load necessary libraries
library(data.table)
library(GenomicAlignments)

# Define the base directory
base_dir <- "/home/hooimin/lu2024-17-19/RAAV-72"

# List all directories in base_dir
directories <- list.dirs(path = base_dir, full.names = TRUE, recursive = FALSE)

# Debugging output: Check the directories found
cat("Directories found:\n")
print(directories)

# Loop through each directory
for (dir in directories) {
  # Extract the base name of the subdirectory
  base_name <- basename(dir)

  # Skip the directory if it is "00_multiqc_result"
  if (base_name == "00_multiqc_result") {
    cat("Skipping directory:", dir, "\n")
    next
  }
  
  # List RDS files matching the pattern in the current directory
  rds_files <- list.files(
    path = file.path(dir, "06_traceBack", "04_p007_found"),
    pattern = "\\.rds$", 
    full.names = TRUE
  )

  # Debugging output: Check the RDS files found
  cat("RDS files found in", dir, ":\n")
  print(rds_files)

  # Initialize an empty list to store GAlignments objects
  alignments_list <- list()
  
  # Read and concatenate the files
  for (file in rds_files) {
    data <- readRDS(file)
    # Debugging output: Check the structure of the data
    cat("Structure of data read from", file, ":\n")
    print(str(data))
    if (inherits(data, "GAlignments")) {
      alignments_list[[length(alignments_list) + 1]] <- data
    } else {
      cat("Skipping invalid data in file:", file, "\n")
    }
  }
  
  # Concatenate all GAlignments objects into one
  if (length(alignments_list) > 0) {
    concatenated_data <- do.call(c, alignments_list)
  } else {
    cat("No valid data found in directory:", dir, "\n")
    next
  }

  # Debugging output: Check the concatenated data
  cat("Concatenated data for", dir, ":\n")
  print(concatenated_data)
  
  # Define the output file path
  output_file <- file.path(dir, "06_traceBack", "04_p007_found", paste0("found.p007_72_", base_name, ".rds"))
  
  # Check if the file already exists and print a message
  if (file.exists(output_file)) {
    cat("File", output_file, "already exists and will be overwritten.\n")
  }
  
  # Save the concatenated data to a new file in the same directory
  saveRDS(concatenated_data, file = output_file)
  
  # Debugging output: Confirm the file has been saved
  cat("Saved concatenated data to", output_file, "\n")
}