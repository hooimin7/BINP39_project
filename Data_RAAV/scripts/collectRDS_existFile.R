#' ---
#' title: "Normalize counts"
#' author: "Hooi Min Tan Grahn"
# Clear the workspace
rm(list = ls())

# Set the working directory to the base directory
setwd("/home/hooimin/lu2024-17-19/RAAV-72")

# Load necessary libraries
library(data.table)
library(fs)  # For file system operations

# Define the base directory
base_dir <- "/home/hooimin/lu2024-17-19/RAAV-72"

# Find all .rds files in the base directory and its subdirectories
rds_files <- list.files(path = base_dir, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)

# Define the output folder path
output_folder <- file.path("../Data_RAAV/p006_data/02_rds")

# Ensure the output folder exists
dir_create(output_folder)

# Copy each .rds file to the output folder
for (file in rds_files) {
  dest_file <- file.path(output_folder, basename(file))
  
  # Check if the file already exists
  if (file_exists(dest_file)) {
    # Generate a new file name to avoid overwriting
    base_name <- tools::file_path_sans_ext(basename(file))
    ext <- tools::file_ext(file)
    counter <- 1
    new_dest_file <- file.path(output_folder, paste0(base_name, "_new_", counter, ".", ext))
    
    # Increment the counter until a non-existing file name is found
    while (file_exists(new_dest_file)) {
      counter <- counter + 1
      new_dest_file <- file.path(output_folder, paste0(base_name, "_new_", counter, ".", ext))
    }
    
    file_copy(file, new_dest_file)
  } else {
    file_copy(file, dest_file)
  }
}