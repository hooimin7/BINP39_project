#' ---
#' title: "Pairwise sample analysis output"
#' author: "Tomas Bjorklund"
#' edited by: "Jaro Steindorff and Hooi Min"
#' output:  
#' ---

#  This will make R-functions such as library() and install.packages() use this directory:
.libPaths(c('~/MyRextensions', .libPaths()))
# Load necessary libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(dplyr))

# Set library paths
.libPaths(c('~/MyRextensions', .libPaths()))

# Start time
strt1 <- Sys.time()

# Read in the all samples data table
all.samples <- readRDS("03_normalize/allSamplesDataTable.RDS")
print(head(all.samples))

# Print the structure of the data table for verification
print("Structure of all.samples:")
str(all.samples)

all.samples <- all.samples %>%
  mutate(across(-c(Group, GeneName), ~ ifelse(is.na(.), "", .)))

# Print the first few rows of the subset data table for verification
print("First few rows of all.samples:")
print(head(all.samples))

# Define the plotting function
plotPair <- function(topSample, bottomSample, size.bin = 1, winWidth = 1, NormalizePlot = TRUE, PlotBC = FALSE) {
  # Create the fill values with correct names
  fill.values <- eval(parse(text = paste("c(", 
                                         paste0(topSample, "= rgb(38,64,135, maxColorValue = 255)"), ", ",
                                         paste0(bottomSample, "= rgb(157,190,217, maxColorValue = 255)"), ")")))
  # Set the key of all samples to the group
  setkey(all.samples, Group)
  # Select the two samples
  select.samples <- all.samples[J(names(fill.values))]

  # Print fill values and group levels for debugging
  print("Fill values:")
  print(names(fill.values))
  print(fill.values)

  print("Group levels in data:")
  print(unique(select.samples$Group))
  print(levels(select.samples$Group))
  
  # Ensure that the Group column contains the correct levels
  if (is.null(unique(select.samples$Group))) {
    stop("Group levels in data are NULL. Ensure that the Group column contains the correct levels.")
  }

  # Adjust the RNA count to log2 and add 1 in case there is a 0
  select.samples[, RNAcount := log2(RNAcount + 1)]

  # Ensure that the Group column contains the correct levels
  print("samples in Group")
  print(unique(select.samples$Group))
  
  # If the window width is greater than 0
  if (winWidth > 0) {
    # Create an ordered table of the samples
    setorder(select.samples, Group, GeneName, start, width)

    windowTable <- select.samples[, c("GeneName", "start", "width"), with = FALSE]
    windowTable <- unique(windowTable, by = c("GeneName", "start", "width"))
    windowTable <- windowTable[, (seq(width - winWidth + 1) + start - 1), by = c("GeneName", "start", "width")]
    setnames(windowTable, "V1", "winStart")
    windowTable[, winEnd := winStart + winWidth - 1]
    setkeyv(windowTable, c("GeneName", "start", "width"))
    setkeyv(select.samples, c("GeneName", "start", "width"))
    select.samples.windowBin <- select.samples[windowTable, allow.cartesian = TRUE]
    select.samples.windowBin[, `:=`(AAproc, winStart/seqlength * 100)]
    
    setkey(select.samples.windowBin, Group)
    select.samples.windowBin <- select.samples.windowBin[J(names(fill.values))] # Select the two compared groups
    setkeyv(select.samples.windowBin, c("Group", "GeneName", "winStart", "winEnd"))

    # Ensure the transformation worked
    print("select.samples.windowBin")
    print(head(select.samples.windowBin))

    # Ensure BC and LUTnrs are character vectors
    if (!is.character(select.samples.windowBin$BC)) {
      stop("BC column is not a character vector")
    }
    if (!is.character(select.samples.windowBin$LUTnrs)) {
      stop("LUTnrs column is not a character vector")
    }

    select.samples.windowBin <- select.samples.windowBin[, list(
      Overlaps = .N,
      seqlength = min(seqlength),
      AAproc = min(AAproc),
      BC = paste(unlist(strsplit(BC, ",")), collapse = ","),  # Concatenate BC strings
      LUTnrs = paste(unlist(strsplit(LUTnrs, ",")), collapse = ","),  # Concatenate LUTnrs strings
      RNAcount = sum(RNAcount)
    ), by = c("Group", "GeneName", "winStart", "winEnd")]
    
    plot.data.dt <- unique(select.samples.windowBin, by = c("Group", "GeneName", "winStart", "winEnd"))

  } else {
    plot.data.dt <- data.table::copy(select.samples)
  }
    
  # Ensure the transformation worked
  print("Transformed select.samples.windowBin")
  print(head(select.samples.windowBin))

  print("plot.data.dt")
  print(head(plot.data.dt))

  # Print intermediate values for debugging
  print("Debugging Values:")
  print(head(plot.data.dt$BC))
  print(head(plot.data.dt$LUTnrs))

  # Binning of data
  FullLength <- 100
  position <- seq(0, FullLength, size.bin)
  plot.data.dt[, bin := findInterval(AAproc, position)]
  
  # Create plot.data.bin with the required calculations
  plot.data.bin <- plot.data.dt[, list(
    .N,  # Number of rows in each group
    seqlength = min(seqlength),  # Minimum sequence length in each group
    AAproc = position[findInterval(mean(AAproc), position)],  # Position of the amino acid
    BCsum = length(table(unlist(strsplit(paste(BC, collapse = ","), ",")))),  # Frequency table of unique barcodes
    LUTnrs = paste(unique(names(table(unlist(strsplit(paste(LUTnrs, collapse = ","), ","))))), collapse = ","),  # Unique LUT numbers
    NormCount = sum(RNAcount) / seqlength * FullLength  # Normalized RNA count
  ), by = c("Group", "GeneName", "bin")]

  print("Intermediate plot.data.bin values:")
  print(head(plot.data.bin))

  plot.data.bin <- unique(plot.data.bin, by = c("Group", "GeneName", "bin"))
  
  print("Final plot.data.bin values:")
  print(head(plot.data.bin))

  # Filtration parameters
  if (NormalizePlot) {
    plot.data.bin[Group == names(fill.values)[1], NormCount := NormCount / max(NormCount)]
    plot.data.bin[Group == names(fill.values)[2], NormCount := NormCount / max(NormCount)]
  }

  # Flip the values for the second group
  plot.data.bin[Group == names(fill.values)[2], NormCount := NormCount * -1]
  
  # Ensure the normalization worked
  print("Normalized plot.data.bin values:")
  print(head(plot.data.bin))

  # Check for NA values in Group column
  if (any(is.na(plot.data.bin$Group))) {
    stop("NA values found in Group column after normalization.")
  }

  # Ensure Group levels are correct
  print("Group levels after normalization:")
  print(unique(plot.data.bin$Group))

  #===================
  # Sort and select top samples
  #===================
  # create a copy of the select samples that we will use to sort and select the top samples
    select.samples.binPos <- select.samples
  # set the key to the group and sequence
    setkeyv(select.samples.binPos,c("Group", "Sequence"))
  # set the order of the samples
    setorder(select.samples.binPos,Group,Sequence,GeneName)
  # select the unique samples. Due to key, this removes replicates if identical sequence mapped to multiple genes
    select.samples.binPos <- unique(select.samples.binPos, by=c("Group","Sequence")) 

    print("unique select bin")
    print(head(select.samples.binPos))
    print("Group levels after unique select bin:")
    print(unique(select.samples.binPos$Group))

  # set the key to the group, category, gene name and AA
    setkeyv(select.samples.binPos,c("Group","Category","GeneName","AA"))

  # Ensure BC and LUTnrs are character vectors
    if (!is.character(select.samples.binPos$BC)) {
      stop("BC column is not a character vector")
    }
    if (!is.character(select.samples.binPos$LUTnrs)) {
      stop("LUTnrs column is not a character vector")
    }

  # Create new columns with the number of BC, the normalized count, the number of animals, the LUT numbers, the main structure, and the number of mismatches
    select.samples.binPos[, `:=`(
      BCcount = length(table(unlist(strsplit(paste(BC, collapse = ","), ",")))),  # Frequency table of unique barcodes
      NormCount = sum(NormCount),
      LUTnrs = paste(unique(names(table(unlist(strsplit(paste(LUTnrs, collapse = ","), ","))))), collapse = ",")  # Unique LUT numbers
    ), by = key(select.samples.binPos)]
  
  # Removes duplicate rows, keeping only unique combinations of values in the columns "Group", "NormCount", and "LUTnrs".
    select.samples.binPos <- unique(select.samples.binPos, by=c("Group","NormCount","LUTnrs"))
  # Selects and retains only the specified columns and removing any columns not that are not listed.
    select.samples.binPos <- select.samples.binPos[,c("Group","GeneName","AA","NormCount",
                                                    "BCcount","LUTnrs"), with = FALSE]
  
  # Print the resulting data table for verification
    print("select.samples.binPos")
    print(head(select.samples.binPos, 10))
    print("Group levels after select.samples.binPos:")
    print(unique(select.samples.binPos$Group))

  # set the key to group
    setkey(select.samples.binPos,Group)
  # select the top 25 samples
    select.samples.top <- select.samples.binPos[, head(.SD, 25), by=Group]

   # Print the resulting data table for verification
    print("select.samples.top")
    print(head(select.samples.top, 10))

    # Extract the gene names of the top 25 samples
    top_genes <- unique(select.samples.top$GeneName)

    # Print the top 25 gene names for verification
    print("Top 25 gene names:")
    print(top_genes)

    # Filter plot.data.bin for the top 25 genes
    plot.data.bin.filtered <- plot.data.bin[GeneName %in% top_genes]

    # Print the filtered plot data for verification
    print("Filtered plot.data.bin values:")
    print(head(plot.data.bin.filtered))

    # Create a list to store the plots
    list_plots <- vector('list', length(top_genes))

    # Create individual plots for each gene
    for (i in seq_along(top_genes)) {
        gene <- top_genes[i]
        plot.data.gene <- plot.data.bin.filtered[GeneName == gene]
        list_plots[[i]] <- ggplot(plot.data.gene, aes(x = AAproc, y = NormCount, fill = Group)) +
            geom_bar(stat = "identity", position = "identity") +
            theme_bw() +
            scale_fill_manual(values = fill.values) +
            labs(title = paste("Gene:", gene), x = "Amino Acid Position", y = "Normalized Count") +
            theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
                  legend.position = "bottom",
                  legend.spacing = unit(0, "cm"),
                  legend.key.height = unit(0, "cm"),
                  plot.background = element_rect(fill = "white"),
                  axis.text = element_text(size = rel(0.45)),
                  axis.ticks = element_line(linewidth = rel(0.5)),
                  axis.ticks.length = unit(0.05, "cm"),
                  strip.text.x = element_text(size = rel(0.5), colour = "black", angle = 0, lineheight = 3, vjust = -20),
                  strip.background = element_blank(),
                  panel.spacing.y = unit(0, "cm"),
                  panel.spacing.x = unit(0, "cm"))
    }

    # Save each plot individually
    for (i in seq_along(list_plots)) {
    ggsave(paste0("04_plot/top25_AAV_lib/plot_", top_genes[i], ".png"), plot = list_plots[[i]], width = 8, height = 6, dpi = 300, units = "in")
    }

    # Optionally, combine all plots into a single PDF
    pdf("04_plot/top25_AAV_lib/combined_plots_AAVlib.pdf", width = 8, height = 6)
    for (i in seq_along(list_plots)) {
        print(list_plots[[i]])
    }
    dev.off()
}

# Generate the plot data using the plotPair function
plot_data <- plotPair("mRNA_All", "library", PlotBC = FALSE)

# Define fill.values using rgb function to generate hexadecimal color values
fill.values <- c(
  "mRNA_All" = rgb(38, 64, 135, maxColorValue = 255),
  "library" = rgb(157, 190, 217, maxColorValue = 255)
)

# To see the warnings, use:
warnings()