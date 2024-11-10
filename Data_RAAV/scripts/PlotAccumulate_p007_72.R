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
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(future.apply))

# Set the working directory to the base directory
setwd("/home/hooimin/lu2024-17-19/Data_RAAV/p007_data/03_normalize")

# Define the base directory
base_dir <- "/home/hooimin/lu2024-17-19/Data_RAAV/p007_data/03_normalize"

# Set the future plan
plan(multicore, workers = 32)

# Start time
strt1 <- Sys.time()

# List all directories in base_dir
directories <- list.dirs(path = base_dir, full.names = TRUE, recursive = FALSE)

# Initialize an empty list to accumulate data
accumulated_data <- list()

# Loop through each directory
for (dir in directories) {
  # List RDS files matching the pattern in the current directory
  rds_files <- list.files(path = dir, pattern = "_normalized_p007_72.*\\.RDS$", full.names = TRUE)
  
  # Process each RDS file
  for (rds_file in rds_files) {
    # Read the RDS file
    all.samples <- readRDS(rds_file)
    
    all.samples <- all.samples %>%
      mutate(across(-c(Group, GeneName), ~ ifelse(is.na(.), "", .)))
    
    # Accumulate data
    accumulated_data[[length(accumulated_data) + 1]] <- all.samples
  }
}

# Combine all accumulated data into a single data.table
combined_data <- rbindlist(accumulated_data)

# Print the combined data for verification
print("Combined data:")
print(head(combined_data))

# Set the key of combined data to the group
setkey(combined_data, Group)

# Select the two samples
select.samples <- combined_data[J(c("mRNA_72", "library"))]

# # Print the first few rows of the selected samples for verification
# print("First few rows of select.samples:")
# print(head(select.samples))

# Ensure that the Group column contains the correct levels
if (is.null(unique(select.samples$Group))) {
  stop("Group levels in data are NULL. Ensure that the Group column contains the correct levels.")
}

# # Adjust the RNA count to log2 and add 1 in case there is a 0
# select.samples[, RNAcount := log2(RNAcount + 1)]

# Define the plotting function
plotPair <- function(topSample, bottomSample, size.bin = 1, winWidth = 1, NormalizePlot = TRUE, PlotBC = FALSE) {
      # Create the fill values with correct names
      fill.values <- eval(parse(text = paste("c(", 
                                             paste0(topSample, "= rgb(38,64,135, maxColorValue = 255)"), ", ",
                                             paste0(bottomSample, "= rgb(157,190,217, maxColorValue = 255)"), ")")))  
  # Set the key of all samples to the group
  setkey(select.samples, Group)
  # Select the two samples
  select.samples <- select.samples[J(c(topSample, bottomSample))]

  # # Print fill values and group levels for debugging
  # print("Fill values:")
  # print(names(fill.values))
  # print(fill.values)

  # print("Group levels in data:")
  # print(unique(select.samples$Group))
  # print(levels(select.samples$Group))
  
  # Ensure that the Group column contains the correct levels
  if (is.null(unique(select.samples$Group))) {
    stop("Group levels in data are NULL. Ensure that the Group column contains the correct levels.")
  }

  # Adjust the RNA count to log2 and add 1 in case there is a 0
  select.samples[, RNAcount := log2(RNAcount + 1)]

  # # Ensure that the Group column contains the correct levels
  # print("samples in Group")
  # print(unique(select.samples$Group))
  
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
    select.samples.windowBin <- select.samples.windowBin[J(c(topSample, bottomSample))] # Select the two compared groups
    setkeyv(select.samples.windowBin, c("Group", "GeneName", "winStart", "winEnd"))

    # # Ensure the transformation worked
    # print("select.samples.windowBin")
    # print(head(select.samples.windowBin))

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

  # print("Intermediate plot.data.bin values:")
  # print(head(plot.data.bin))

  plot.data.bin <- unique(plot.data.bin, by = c("Group", "GeneName", "bin"))
  
  # print("Final plot.data.bin values:")
  # print(head(plot.data.bin))

  # Filtration parameters
  if (NormalizePlot) {
    plot.data.bin[Group == topSample, NormCount := NormCount / max(NormCount)]
    plot.data.bin[Group == bottomSample, NormCount := NormCount / max(NormCount)]
  }

  # Flip the values for the second group
  plot.data.bin[Group == bottomSample, NormCount := NormCount * -1]
  
  # # Ensure the normalization worked
  # print("Normalized plot.data.bin values:")
  # print(head(plot.data.bin))

  # Check for NA values in Group column
  if (any(is.na(plot.data.bin$Group))) {
    stop("NA values found in Group column after normalization.")
  }

  # # Ensure Group levels are correct
  # print("Group levels after normalization:")
  # print(unique(plot.data.bin$Group))

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
  

  # set the key to group
    setkey(select.samples.binPos,Group)
  # select the top 100 samples
    select.samples.top <- select.samples.binPos[, head(.SD, 100), by=Group]

    # Extract and process topSample and bottomSample
    topSample <- select.samples.top[J(topSample)]
    bottomSample <- select.samples.top[J(bottomSample)]

    # Remove the group column and set the gene name as the column name
    topSample[, c("Group") := NULL]
    bottomSample[, c("Group") := NULL]
    
    # Extract the gene names of the top 100 samples
    top_genes <- unique(select.samples.top$GeneName)

    # Print the top 100 gene names for verification
    print("Top 100 gene names:")
    print(top_genes)

    # Filter plot.data.bin for the top 100 genes
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
    ggsave(paste0("../04_plot/top25_AAV_lib_test/combined_plot_72AAV_", top_genes[i], ".png"), plot = list_plots[[i]], width = 8, height = 6, dpi = 300, units = "in")
    }

    # Optionally, combine all plots into a single PDF
    pdf("../04_plot/top25_AAV_lib_test/combined_plots_72AAV.pdf", width = 8, height = 6)
    for (i in seq_along(list_plots)) {
        print(list_plots[[i]])
    }
    dev.off()
}

# Call the plotPair function with appropriate arguments
plotPair("mRNA_72", "library")

# Define fill.values using rgb function to generate hexadecimal color values
fill.values <- c(
  "mRNA_72" = rgb(38, 64, 135, maxColorValue = 255),
  "library" = rgb(157, 190, 217, maxColorValue = 255)
)

# To see the warnings, use:
warnings()

# Print total analysis time
print("Total analysis time:")
print(Sys.time() - strt1)

# Print session info
devtools::session_info()