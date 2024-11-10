#' ---
#' title: "Normalize counts"
#' author: "Hooi Min Tan Grahn"
# Clear the workspace
rm(list = ls())

# Set the working directory to the base directory
setwd("/home/hooimin/lu2024-17-19/RAAV-60")

# Load necessary libraries
library(data.table)
library(fs)  # For file system operations

# Define the base directory
base_dir <- "/home/hooimin/lu2024-17-19/RAAV-60"

# Find all .rds files in the base directory and its subdirectories
rds_files <- list.files(path = base_dir, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)

# Define the output folder path
output_folder <- file.path("../Data_RAAV/p006_data/02_rds")

# Copy each .rds file to the output folder
for (file in rds_files) {
  file_copy(file, file.path(output_folder, basename(file)))
}
