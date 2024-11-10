#' ---
#' title: "Normalize counts"
#' author: "Hooi Min Tan Grahn"
#' ---

# Clear the workspace
rm(list = ls())

suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(future.apply))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(BiocGenerics))
suppressPackageStartupMessages(library(GenomicRanges)) 

#' Generate load list and grouping names
#' ============================
strt <- Sys.time()

# Base directory
base_dir <- "/home/hooimin/lu2024-17-19/Data_RAAV/p006_data"


# Load the list of files to normalize
in.names.all <- list.files(path = base_dir, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)

# Read the combined log list into a data frame
load.list <- fread("01_input/cleanloglist.txt", header = FALSE, sep = "\t", 
                   stringsAsFactors = FALSE, fill = TRUE)


load.list <- load.list[-1, ] # Remove the first row

# Assign column names
setnames(load.list, c("Name", "GroupName"))


# print("load.list after set names")
# print(load.list)

# Select the files to load
select.Cases <- c(unlist(sapply(load.list$Name, function(x) grep(x, in.names.all),
                                simplify = TRUE)))

select.Cases <- unique(select.Cases)

# Load the files
(in.names.all <- in.names.all[select.Cases])

# Ensure select.Cases has names
names(select.Cases) <- load.list$Name[select.Cases]

# Create the grouping data frame
grouping <- data.frame(
  Sample = gsub("-", "_", gsub("found.|(02_rds/)|(.rds)", "", basename(in.names.all))),
  Group = load.list$GroupName[match(names(select.Cases), load.list$Name)],
  stringsAsFactors = FALSE
)

print("grouping")
print(head(grouping))

# Function to load the RDS files and handle different formats
loadRDS <- function(in.name) {
  this.sample <- readRDS(in.name)
  
  # Print the class of this.sample for debugging
  print(paste("Loading:", in.name))
  print(paste("Class of this.sample:", class(this.sample)))
  
  this.name <- gsub("-", "_", gsub("found.|(02_rds/)|(.rds)", "", basename(in.name)))
  this.group <- grouping[match(this.name, grouping$Sample), "Group"]
  
  if (inherits(this.sample, "GAlignments")) {
    # Handle GAlignments object
    mcols(this.sample) <- cbind(mcols(this.sample),
                                data.frame(Sample = this.name, Group = this.group,
                                           stringsAsFactors = FALSE))
  } else if (inherits(this.sample, c("data.table", "data.frame"))) {
    # Handle data.table or data.frame object
    this.sample[, Sample := this.name]
    this.sample[, Group := this.group]
  } else {
    warning(paste("Skipping file:", in.name, "- Unsupported object type"))
    return(NULL)
  }
  
  return(this.sample)
}

# # Load all RDS files specified in 'in.names.all' into a list of R objects, filtering out NULL results
all.samples <- lapply(in.names.all, loadRDS)
all.samples <- Filter(Negate(is.null), all.samples)

print("all.samples after loading")
print(head(all.samples))

# Combine the list of R objects into a single GAlignmentsList object. Using unlist() to flatten the list of objects before combining
all.samples <- GAlignmentsList(unlist(all.samples))

print("all.samples after GA")
print(head(all.samples))

# Rename the samples into unique and R valid names
names(all.samples) <- make.names(names(all.samples), unique = TRUE)

print("after unique names")
print(all.samples)

# create a data.table with the seqnames and seq length and set the key to seqnames  
length.Table <- data.table(seqnames = names(seqlengths(all.samples)), seqlength = seqlengths(all.samples),
key = "seqnames")

print(length.Table)

# create a data.table from all samples and set the key to seqnames    
all.samples <- data.table(as.data.frame(all.samples), key = "seqnames")

# Inspect the structure of the original all.samples data frame
str(all.samples)

print(" seqnames")
print(all.samples)

# remove unnecessary columns
all.samples[, `:=`(c("strand", "qwidth", "cigar", "njunc", "end", "Frequency", "i.Name"), NULL)]

print("remove 5 col")
print(all.samples)

# merge all samples with the length table to get the sequence length of every sequence
all.samples <- all.samples[length.Table] #A data.table merge to match seqlengths to their respective seqnames

print("after merge with length")
print(all.samples)

# create new columns and assign the values based on the seqnames that was split into its parts
all.samples[, `:=`(c("Category", "Protein", "Origin", "Extra", "Number", "GeneName"),
tstrsplit(seqnames, ",", fixed = TRUE))]

print("after split")
print(all.samples)

# Remove unnecessary columns
all.samples[, `:=`(c("seqnames", "Protein", "Origin", "Extra", "Number"), NULL)]

# Change the GeneName column to hold the actual gene name
all.samples[, `:=`(GeneName, gsub("#\\d+", "", GeneName))]

print("Change genename")
print(all.samples)


#' Normalizing read counts to correct for variable read depth
#' ============================
# setting the key to the group
setkey(all.samples,Group)

# get rid of single read samples
all.samples <- all.samples[RNAcount>1,]

print("remove single read")
print(all.samples)

# store the read counts grouped by the group
readCounts <- all.samples[,list(GroupCount=sum(RNAcount)), by="Group"]

print("stored read count")
print(readCounts)

# Normalize the group reads by the max Group reads set values between 0 and 1 where 1 is the max value
readCounts[,GroupCount:=GroupCount/max(GroupCount)]

print("Normalized read count")
print(readCounts)

# Extract only the definitive Barcodes which have at least 2 reads per pair of barcode & LUTnr
all.samples <- all.samples[Mode=="Def"]

print("samples Def")
print(all.samples)

# set the key to Group
setkey(readCounts,Group)
# Add the normalized group count to the all samples data table
all.samples <- all.samples[readCounts]

print("normalized group count")
print(all.samples)

# normalize the read counts by the Group Count to correct for variable read depth
all.samples[,RNAcount:=RNAcount/GroupCount]

print("normalized read depth")
print(all.samples)

# set the key to the Group
setkey(all.samples,Group)

# get all the AAV samples by removing the DNA samples
total.AAV.samples <- all.samples[Group!="library"]
# print the number of AAV samples
print(nrow(total.AAV.samples))
print("Remove DNA")
print(head(total.AAV.samples))

# get the transported AAV samples
transported.AAV.samples <- total.AAV.samples[grepl("AAV",total.AAV.samples$Group)]
# print the number of AAV samples
print(nrow(transported.AAV.samples))
print("AAV")
print(head(transported.AAV.samples))

# get the transported 53 tissue AAV samples
transported.AAV.samples.53 <- total.AAV.samples[grepl("RAAV-53",total.AAV.samples$Group)]
# print the number of 53 tissue AAV samples
print(nrow(transported.AAV.samples.53))
# set the group name to mRNA_53
transported.AAV.samples.53[,Group := "mRNA_53"]
print("mRNA_53")
print(head(transported.AAV.samples.53))

# get the transported 72 tissue AAV samples
transported.AAV.samples.72 <- total.AAV.samples[grepl("RAAV-72",total.AAV.samples$Group)]
# print the number of 72 tissue AAV samples
print(nrow(transported.AAV.samples.72))
# set the group name to mRNA_72
transported.AAV.samples.72[,Group := "mRNA_72"]
print("mRNA_72")
print(head(transported.AAV.samples.72))

# set the group name to mRNA_All
total.AAV.samples[,Group := "mRNA_All"]

# recombine all samples with their new group names
all.samples <- rbind(all.samples,total.AAV.samples,transported.AAV.samples.53,transported.AAV.samples.72)

print("after rbind")
print(all.samples)

# remove the data.tables that were create before
rm(total.AAV.samples,transported.AAV.samples.53,transported.AAV.samples.72)

# set the key of all samples to a vector of the columns
setkeyv(all.samples,c("Group","Category","GeneName","start","width","Sequence","seqlength"))

# now combine all information about identical fragments that are found with different barcodes
all.samples <- all.samples[,j=list(mCount=sum(mCount),
                                   tCount=sum(tCount),
                                   BC=paste(unique(BC), collapse = ","),
                                   Animals=paste(unique(Sample), collapse = ","),
                                   LUTnrs=paste(unique(LUTnr), collapse = ","),
                                   RNAcount=sum(RNAcount),
                                   NormCount=log2(sum(RNAcount)+1)*.N),
                           by=c("Group","Category","GeneName","start","width","Sequence","seqlength")]

print("Identical frag")
print(all.samples)

# adjust the start,width and seqlength from DNA based to AA based
all.samples[,start:=floor((start+2)/3)]
print("start")
print(all.samples)

all.samples[,width:=ceiling((width)/3)]
print("width")
print(all.samples)

all.samples[,seqlength:=ceiling(seqlength/3)]
print("seqlength")
print(all.samples)

# calculate the absolute AA position of the middle of the fragment 
all.samples[,AA:=floor(start+(width/2))]
print("AA pos")
print(all.samples)

# calculate the relative AA position of the middle of the fragment in the whole sequence
all.samples[,AAproc:=AA/seqlength*100]
print("relative AA pos")
print(all.samples)

#Change the default behavior to induce start codons and Methionine
GENETIC_CODE_ALT <- GENETIC_CODE
attr(GENETIC_CODE_ALT, "alt_init_codons") <- c("TAA","TAG")

# Define the translation function
translate_function <- function(sequence) {
  as.character(Biostrings::translate(DNAString(sequence), genetic.code = GENETIC_CODE_ALT, if.fuzzy.codon = "solve"))}

# Set up parallel processing
available_cores <- availableCores()  # Get the number of available cores
plan(multisession, workers = available_cores)  # Adjust number of workers based on available cores

# Apply translation function in parallel
all.samples[, Peptide := future_lapply(Sequence, translate_function)]

# change the Peptide column to character
all.samples[, Peptide := as.character(Peptide)]

# Optionally, reset the parallel plan
plan(sequential)

# print the number of samples in the all samples data table
print(nrow(all.samples))

# save the whole data table as all Samples DATA table
saveRDS(all.samples, file="03_normalize/allSamplesDataTable.RDS")

print("Total execution time:")
print(Sys.time()-strt)

