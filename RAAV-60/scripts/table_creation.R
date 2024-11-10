rm(list=ls()) # Clear workspace

# Load necessary package
library(knitr)

# Read the CSV file
csv_file_path <- 'sorted_fragment_align_7_sub.csv'
data <- read.csv(csv_file_path)

# Create a Markdown table
markdown_table <- kable(data, format = "markdown")

# Define the output file path
output_file_path <- 'table_7.md'

# Save the Markdown table to a file
writeLines(markdown_table, output_file_path)
