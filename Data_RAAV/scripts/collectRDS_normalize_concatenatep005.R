# Clear the workspace
rm(list = ls())

# Set the working directory to the base directory
setwd("/home/hooimin/lu2024-17-19/Data_RAAV/p005_data/03_normalize")

# Load necessary libraries
library(data.table)

# Define the base directory
base_dir <- "/home/hooimin/lu2024-17-19/Data_RAAV/p005_data/03_normalize"

# List all directories in base_dir
directories <- list.dirs(path = base_dir, full.names = TRUE, recursive = FALSE)

# Debugging output: Check the directories found
cat("Directories found:\n")
print(directories)

# Loop through each directory
for (dir in directories) {
  # Extract the base name of the subdirectory
  base_name <- basename(dir)
  
  # List RDS files matching the pattern in the current directory
  # rds_files <- list.files(path = dir, pattern = "allSamplesDataTable_p005_53_[0-9]+_chunk[0-9]+\\.RDS$", full.names = TRUE)
  # rds_files <- list.files(path = dir, pattern = "allSamplesDataTable_p005_72_S[0-9]+_chunk[0-9]+\\.RDS$", full.names = TRUE)
  rds_files <- list.files(
  path = dir, 
  pattern = "allSamplesDataTable_p005_(53_[0-9]+|72_S[0-9]+)_chunk[0-9]+\\.RDS$", 
  full.names = TRUE
  )

  # Initialize an empty list to store data frames
  data_list <- list()
  
  # Read and concatenate the files
  for (file in rds_files) {
    data <- readRDS(file)
    data_list[[length(data_list) + 1]] <- data
  }
  
  # Concatenate all data frames into one
  concatenated_data <- rbindlist(data_list)

  print("processing")
  
  # Define the output file path
  output_file <- file.path(dir, paste0("allSamplesDataTable_p005_", base_name, ".RDS"))
  
  # Check if the file already exists and print a message
  if (file.exists(output_file)) {
    cat("File", output_file, "already exists and will be overwritten.\n")
  }
  
  # Save the concatenated data to a new file in the same directory
  saveRDS(concatenated_data, file = output_file)
  
}