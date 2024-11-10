#' ---
#' title: "Pairwise sample analysis output"
#' author: "Tomas Bjorklund"
#' edited by: "Hooi Min"
#' output: html_document
#' ---

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(devtools))

# Start time
strt1 <- Sys.time()

# Define the base directory
base_dir <- "/home/hooimin/lu2024-17-19/Data_RAAV/p005_data"

LUT.dna <- fread("../../RAAV-60/p005/04_blast/LUTdna.csv")

# Define the path to the 03_normalize subdirectory
normalize_dir <- file.path(base_dir, "03_normalize")

# List all subdirectories within the 03_normalize directory
subdirectories <- list.dirs(path = normalize_dir, full.names = TRUE, recursive = TRUE)

# Filter out the 03_normalize directory itself
subdirectories <- subdirectories[subdirectories != normalize_dir]

# Print the result
print(subdirectories)

# Define patterns for the two sets of files
patterns <- c(
  "allSamplesDataTable_p005_(26_S3|27_S4|28_S5|30_S6|31_S7|32_S8|34_S9|35_S10|36_S11|37_S12|38_S13|39_S14|42_S15|43_S16|44_S17|46_S18|49_S19|50_S20|51_S21|52_S22)\\.RDS$",
  "allSamplesDataTable_p005_(S26_S1|S28_S2|S54_S3|S55_S4|S56_S5|S57_S6|S58_S7|S59_S8)\\.RDS$"
)

# Loop through each subdirectory and process files
for (subdir in subdirectories) {
  base_name <- basename(subdir)
  
  for (pattern in patterns) {
    data <- process_rds_files(subdir, pattern)
    if (!is.null(data)) {
      save_extracted_sequences(data, subdir, base_name)
    }
  }
}

# Combine all seq.lib_split, seq.AAV, seq.53, and seq.72 from all subdirectories
all_seq.lib_split <- unique(unlist(lapply(subdirectories, function(subdir) {
  readRDS(file.path(subdir, paste0("seq_lib_split_", basename(subdir), ".RDS")))
})))

all_seq.AAV <- unique(unlist(lapply(subdirectories, function(subdir) {
  readRDS(file.path(subdir, paste0("seq_AAV_", basename(subdir), ".RDS")))
})))

all_seq.53 <- unique(unlist(lapply(subdirectories, function(subdir) {
  readRDS(file.path(subdir, paste0("seq_53_", basename(subdir), ".RDS")))
})))

all_seq.72 <- unique(unlist(lapply(subdirectories, function(subdir) {
  readRDS(file.path(subdir, paste0("seq_72_", basename(subdir), ".RDS")))
})))

# Continue with the script
venn.area1 <- length(LUT.dna$LUTnr)
print(venn.area1)
venn.area2 <- length(all_seq.lib_split)
print(venn.area2)
venn.area3 <- length(all_seq.AAV)
print(venn.area3)
venn.area4 <- length(all_seq.53)
print(venn.area4)
venn.area5 <- length(all_seq.72)
print(venn.area5)

isect.53_72 <- length(intersect(all_seq.53, all_seq.72))
print("isect.53_72")
print(isect.53_72)
venn.n12 <- length(intersect(LUT.dna$LUTnr, all_seq.lib_split))
print("venn.n12")
print(venn.n12)
venn.n23 <- length(intersect(all_seq.lib_split, all_seq.AAV))
print("venn.n23")
print(venn.n23)
venn.n13 <- length(intersect(LUT.dna$LUTnr, all_seq.AAV))
print("venn.n13")
print(venn.n13)
venn.n123 <- length(intersect(intersect(LUT.dna$LUTnr, all_seq.lib_split), all_seq.AAV))
print("venn.n123")
print(venn.n123)
isect.AAV_53 <- length(intersect(all_seq.AAV, all_seq.53))
print("isect.AAV_53")
print(isect.AAV_53)
isect.AAV_72 <- length(intersect(all_seq.AAV, all_seq.72))
print("isect.AAV_72")
print(isect.AAV_72)

# Identify LUTnr values in seq.lib that are not in seq.arry
lib_not_in_arry <- setdiff(all_seq.lib_split, LUT.dna$LUTnr)
print("LUTnr values in seq.lib but not in seq.arry")
print(lib_not_in_arry)

output.table <- data.frame(
  NameArray = character(),
  NameLib = character(),
  NameAAV = character(),
  Name53 = character(),
  Name72 = character(),
  ArrayStart = numeric(),
  ArrayEnd = numeric(),
  LibStart = numeric(),
  LibEnd = numeric(),
  AAVStart = numeric(),
  AAVEnd = numeric(),
  tissue53Start = numeric(),
  tissue53End = numeric(),
  tissue72Start = numeric(),
  tissue72End = numeric(),
  stringsAsFactors = FALSE
)
output.table[1:2, 1] <- c("Array", "None")
output.table[1:2, 2] <- c("library", "None")
output.table[1:2, 3] <- c("AAV", "None")
output.table[1:2, 4] <- c("mRNA_53", "None")
output.table[1:2, 5] <- c("mRNA_72", "None")
output.table[1:2, 6:12] <- 0

print("output.table empty")
print(output.table)

output.table$ArrayEnd[1] <- output.table$LibEnd[2] <- output.table$AAVEnd[2] <- 
  output.table$tissue53End[2] <- output.table$tissue72End[2] <- length(LUT.dna$LUTnr)

output.table$LibStart[2] <- output.table$LibEnd[1] <- length(intersect(LUT.dna$LUTnr, all_seq.lib_split))

output.table$AAVStart[2] <- output.table$AAVEnd[1] <- length(intersect(all_seq.lib_split, all_seq.AAV))

output.table$tissue53Start[2] <- output.table$tissue53End[1] <- length(intersect(all_seq.AAV, all_seq.53)) 

output.table$tissue72Start[2] <- output.table$tissue72End[1] <- length(intersect(all_seq.AAV, all_seq.72)) 

print("output.table after filled")
print(output.table)

fill.values <- c(
  Array = rgb(193, 210, 234, maxColorValue = 255),
  library = rgb(117, 160, 207, maxColorValue = 255),
  AAV = rgb(157, 190, 217, maxColorValue = 255),
  mRNA_72 = rgb(117, 160, 207, maxColorValue = 255),
  mRNA_53 = rgb(38, 64, 135, maxColorValue = 255),
  None = rgb(255, 255, 255, maxColorValue = 255, alpha = 0)
)

plot <- ggplot(output.table) + 
  scale_x_continuous(limit = c(0, 10), breaks = c(seq(1, 8, 1)), expand = c(0, 0)) + 
  scale_y_continuous(breaks = c(seq(0, 150000, 15000))) + 
  scale_fill_manual(name = "Library", values = fill.values) + 
  theme(aspect.ratio = 1) +
  geom_rect(data = output.table, aes(fill = Name53, ymax = tissue53End, ymin = tissue53Start, xmax = 7.95, xmin = 6.05)) + 
  geom_rect(data = output.table, aes(fill = Name72, ymax = tissue72End, ymin = tissue53End, xmax = 7.95, xmin = 6.05)) +
  geom_rect(data = output.table, aes(fill = NameAAV, ymax = AAVEnd, ymin = AAVStart, xmax = 5.95, xmin = 4.05)) + 
  geom_rect(data = output.table, aes(fill = NameLib, ymax = LibEnd, ymin = LibStart, xmax = 3.95, xmin = 2.05)) +
  geom_rect(data = output.table, aes(fill = NameArray, ymax = ArrayEnd, ymin = ArrayStart, xmax = 1.95, xmin = 0)) +
  coord_polar(theta = "y")

# Save the plot
ggsave("05_vennplot/circular_venn_plot.png", plot = plot, width = 10, height = 10, units = "in")