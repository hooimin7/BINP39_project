#' ---
#' title: "Normalize counts"
#' author: "Hooi Min Tan Grahn"
# Clear the workspace
rm(list = ls())

library(data.table)

# Set the base directory
base_dir <- "/home/hooimin/lu2024-17-19/Data_RAAV/p006_data/02_rds"

# Set the working directory to the base directory
setwd(base_dir)

# Find all .rds files in the base directory and its subdirectories
rds_files <- list.files(path = base_dir, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)

# Define the output file path
output_file <- file.path("../01_input/cleanloglist.txt")

# Clean the names
cleaned_names <- gsub("-", "_", gsub("found.|(.rds)", "", basename(rds_files)))

# Create a data.table with two columns: Name and GroupName
cleaned_data <- data.table(Name = cleaned_names, GroupName = "")

# Define the group name patterns and assign group names
cleaned_data[grepl("multipleContfragmentsCompleteNoLUT", Name), GroupName := "library"]
cleaned_data[grepl("_53_", Name), GroupName := "RAAV-53"]
cleaned_data[grepl("_72_", Name), GroupName := "RAAV-72"]
cleaned_data[grepl("_AAV", Name), GroupName := "AAV"]

# Save the cleaned data as a TXT file
fwrite(cleaned_data, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)


