#' ---
#' title: "Pairwise sample analysis output"
#' author: "Tomas Bjorklund"
#' edited by: "Hooi Min"
#' output: html_document
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

LUT.dna <- fread("../../RAAV-60/p007/04_blast/LUTdna.csv")

complete.library <- readRDS("03_normalize/allSamplesDataTable_p007.RDS")

setkey(complete.library, Group)

seq.arry <- LUT.dna$LUTnr

print("seq.arry")
print(head(seq.arry))

# Extract unique LUTnr values from the library group
seq.lib <- unique(complete.library[J("library")]$LUTnr)

# Split concatenated LUTnr values in seq.lib into individual values
seq.lib_split <- unlist(strsplit(seq.lib, ","))

print("seq.lib")
print(head(seq.lib))

seq.AAV <- unique(complete.library[J("mRNA_All")]$LUTnr)

print("seq.AAV")
print(head(seq.AAV))

seq.53 <- unique(complete.library[grep("mRNA_53", Group)]$LUTnr)

print("seq.53")
print(head(seq.53))

seq.72 <- unique(complete.library[grep("mRNA_72", Group)]$LUTnr)

print("seq.72")
print(head(seq.72))

venn.area1 <- length(seq.arry)
print(venn.area1)
venn.area2 <- length(seq.lib)
print(venn.area2)
venn.area3 <- length(seq.AAV)
print(venn.area3)
venn.area4 <- length(seq.53)
print(venn.area4)
venn.area5 <- length(seq.72)
print(venn.area5)

isect.53_72 <- length(intersect(seq.53, seq.72))
print("isect.53_72")
print(isect.53_72)
venn.n12 <- length(intersect(seq.arry, seq.lib_split))
print("venn.n12")
print(venn.n12)
venn.n23 <- length(intersect(seq.lib, seq.AAV))
print("venn.n23")
print(venn.n23)
venn.n13 <- length(intersect(seq.arry, seq.AAV))
print("venn.n13")
print(venn.n13)
venn.n123 <- length(intersect(intersect(seq.arry, seq.lib), seq.AAV))
print("venn.n123")
print(venn.n123)
isect.AAV_53 <- length(intersect(seq.AAV, seq.53))
print("isect.AAV_53")
print(isect.AAV_53)
isect.AAV_72 <- length(intersect(seq.AAV, seq.72))
print("isect.AAV_72")
print(isect.AAV_72)

# Identify LUTnr values in seq.lib that are not in seq.arry
lib_not_in_arry <- setdiff(seq.lib_split, seq.arry)
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
output.table$tissue53End[2] <- output.table$tissue72End[2] <- length(seq.arry)

output.table$LibStart[2] <- output.table$LibEnd[1] <- length(intersect(seq.arry, seq.lib_split))

output.table$AAVStart[2] <- output.table$AAVEnd[1] <- length(intersect(seq.lib, seq.AAV))

output.table$tissue53Start[2] <- output.table$tissue53End[1] <- length(intersect(seq.AAV, seq.53)) 

output.table$tissue72Start[2] <- output.table$tissue72End[1] <- length(intersect(seq.AAV, seq.72)) 

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
ggsave("05_vennplot/circular_venn_plot.png", plot = plot, width = 10, height = 10, unit = "in")

# Save the plot as EPS
ggsave("05_vennplot/circular_venn_plot.eps", plot = plot, width = 10, height = 10, unit = "in")