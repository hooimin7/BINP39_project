#' ---
#' title: "Combine several data"
#' author: "Hooi Min Tan Grahn"
# Clear the workspace
rm(list = ls())

# Load necessary libraries
library(data.table)  # Optional

# Base directory
base_dir <- "/home/hooimin/lu2024-17-19/Data_RAAV"

# Read the first file into a data frame
df_RAAV_30 <- fread(file.path(base_dir, "../RAAV-30/p005_AAV_01/06_traceBack/01_table/loglist_30.txt"), sep = "\t", header = TRUE)

# Read the second file into a data frame
df_RAAV_53 <- fread(file.path(base_dir, "../RAAV-53/26_S3_L001_R1_001/06_traceBack/01_table/loglist_53.txt"), sep = "\t", header = TRUE)

df_RAAV_72 <- fread(file.path(base_dir, "../RAAV-72/S26_S1_L001_R1_001/06_traceBack/01_table/loglist_72.txt"), sep = "\t", header = TRUE)

# Concatenate the two data frames
combined_df <- rbind(df_RAAV_30, df_RAAV_53, df_RAAV_72)

# Create a new data frame with the values "p005", "p006", and "p007"
new_rows <- data.frame(Name = c("p005", "p006", "p007"))

# Concatenate the new rows with the combined data frame
combined_df <- rbind(combined_df, new_rows)

# Define the output file name and path
output_dir <- file.path(base_dir, "01_input")
output_file <- file.path(output_dir, "combined_loglist.txt")

# Ensure the output directory exists
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Save the combined_df as a TXT file
fwrite(combined_df, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)


