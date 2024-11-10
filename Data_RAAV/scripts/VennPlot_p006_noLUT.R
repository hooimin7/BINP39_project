#' ---
#' title: "Pairwise sample analysis output"
#' author: "Tomas Bjorklund"
#' edited by: "Hooi Min"
#' output:  
#' ---

# This will make R-functions such as library() and install.packages() use this directory:
.libPaths(c('~/MyRextensions', .libPaths()))

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

LUT.lib <- readRDS("../../RAAV-60/p006/04_blast/02_analysis/multipleContfragmentsCompleteNoLUT.rds")

print(head(LUT.lib))

# Add the Group column with the value "library"
LUT.lib[, Group := "library"]

# Set the key to the Group column
setkey(LUT.lib, Group)

# Extract the Reads column for the "library" group
seq.lib <- unique(LUT.lib[J("library")]$Reads)

print("seq.lib")
print(head(seq.lib))

normalize.library <- readRDS("03_normalize/allSamplesDataTable.RDS")

print(head(normalize.library))

# Set the key to the Group column
setkey(normalize.library, Group)

# Extract the Reads column for the "AAV" group
seq.AAV <- unique(normalize.library[J("AAV")]$Reads)

print("seq.AAV")
print(head(seq.AAV))

# Extract the Reads column for the "mRNA_53" group
seq.53 <- unique(normalize.library[grep("mRNA_53", Group)]$Reads)

print("seq.53")
print(head(seq.53))

# Extract the Reads column for the "mRNA_72" group
seq.72 <- unique(normalize.library[grep("mRNA_72", Group)]$Reads)

print("seq.72")
print(head(seq.72))

# Calculate the lengths of the sequences
venn.area2 <- length(seq.lib)
print(venn.area2)
venn.area3 <- length(seq.AAV)
print(venn.area3)
venn.area4 <- length(seq.53)
print(venn.area4)
venn.area5 <- length(seq.72)
print(venn.area5)

# Calculate intersections
isect.53_72 <- length(intersect(seq.53, seq.72))
print("isect.53_72")
print(isect.53_72)
isect.AAV.53 <- length(intersect(seq.AAV, seq.53))
print("isect.AAV_53")
print(isect.AAV.53)
isect.AAV.72 <- length(intersect(seq.AAV, seq.72))
print("isect.AAV_72")
print(isect.AAV.72)
venn.n23 <- length(intersect(seq.lib, seq.AAV))
print("venn.n23")
print(venn.n23)

# Create the output table
output.table <- data.frame(
  NameLib = character(),
  NameAAV = character(),
  Name53 = character(),
  Name72 = character(),
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
output.table[1:2, 1] <- c("library", "None")
output.table[1:2, 2] <- c("AAV", "None")
output.table[1:2, 3] <- c("mRNA_53", "None")
output.table[1:2, 4] <- c("mRNA_72", "None")
output.table[1:2, 5:12] <- 0

print("output.table empty")
print(output.table)

# Fill the output table with calculated values
output.table$LibEnd[1] <- output.table$AAVEnd[2] <- 
output.table$tissue53End[2] <- output.table$tissue72End[2] <- length(seq.lib)

output.table$AAVStart[2] <- output.table$AAVEnd[1] <- 
length(intersect(seq.lib, seq.AAV))

output.table$tissue53Start[2] <- output.table$tissue53End[1] <- length(intersect(seq.AAV, seq.53)) 

# output.table$tissue53Start[2] <- output.table$tissue53End[1] <- length(intersect(seq.AAV, seq.53)) - length(intersect(seq.AAV, seq.72))

output.table$tissue72Start[2] <- output.table$tissue72End[1] <- length(intersect(seq.AAV, seq.72))

print("output.table after filled")
print(output.table)

# Define fill values for the plot
fill.values <- c(
  library = rgb(193, 210, 234, maxColorValue = 255),
  AAV = rgb(157, 190, 217, maxColorValue = 255),
  mRNA_72 = rgb(117, 160, 207, maxColorValue = 255),
  mRNA_53 = rgb(38, 64, 135, maxColorValue = 255),
  None = rgb(255, 255, 255, maxColorValue = 255, alpha = 0)
)

# Create the plot
plot <- ggplot(output.table) + 
  scale_x_continuous(limit = c(0, 6), breaks = c(seq(1, 6, 1)), expand = c(0, 0)) + 
  scale_y_continuous(breaks = c(seq(0, 100000, 20000))) + 
  scale_fill_manual(name = "Library", values = fill.values) + 
  theme(aspect.ratio = 1) +
  geom_rect(data = output.table, aes(fill = Name53, ymax = tissue53End, ymin = tissue72End, xmax = 5.95, xmin = 4.05)) + 
  geom_rect(data = output.table, aes(fill = Name72, ymax = tissue72End, ymin = tissue72Start, xmax = 5.95, xmin = 4.05)) +
  geom_rect(data = output.table, aes(fill = NameAAV, ymax = AAVEnd, ymin = AAVStart, xmax = 3.95, xmin = 2.05)) + 
  geom_rect(data = output.table, aes(fill = NameLib, ymax = LibEnd, ymin = LibStart, xmax = 1.95, xmin = 0)) +
  coord_polar(theta = "y")

# Save the plot
ggsave("05_vennplot/circular_venn_plot.png", plot = plot, width = 10, height = 10, unit = "in")

# Save the plot as EPS
ggsave("05_vennplot/circular_venn_plot.eps", plot = plot, width = 10, height = 10, unit = "in")