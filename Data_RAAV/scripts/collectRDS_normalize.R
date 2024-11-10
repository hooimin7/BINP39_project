# Clear the workspace
rm(list = ls())

# Set the working directory to the base directory
setwd("/home/hooimin/lu2024-17-19/Data_RAAV/p007_data/03_normalize")

# Load necessary libraries
library(data.table)

# Define the base directory
base_dir <- "/home/hooimin/lu2024-17-19/Data_RAAV/p007_data/03_normalize"

# List files matching the pattern _normalized_xxxx.RDS
rds_files <- list.files(path = base_dir, pattern = "_normalized_.*\\.RDS$", recursive = TRUE, full.names = TRUE)

# Initialize an empty list to store data frames
data_list <- list()

# Read and concatenate the files
for (file in rds_files) {
  data <- readRDS(file)
  data_list[[length(data_list) + 1]] <- data
}

# Concatenate all data frames into one
concatenated_data <- rbindlist(data_list)

# Save the concatenated data to a new file
saveRDS(concatenated_data, file = file.path(base_dir, "allSamplesDataTable_p007.RDS"))

# Print a message indicating success
cat("Concatenated data saved to allSamplesDataTable_p007.RDS\n")