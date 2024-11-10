
#' ---
#' title: "Clustering Polypeptide motifs using the Hammock hidden Markov model peptide clustering"
#' author: "Tomas Bjorklund"
#' edited by: "Hooi Min Tan Grahn"
#' header-includes: \usepackage{longtable}
#' output: 
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.6in
#' fontsize: 9.5pt
#' ---

# Clear the workspace
rm(list = ls())

setwd("/home/hooimin/lu2024-17-19/Data_RAAV/p006_data")

#' This script clusters Polypeptide motifs using the Hammock hidden Markov model peptide clustering.  
suppressPackageStartupMessages(library(knitr))
#+ setup, include=FALSE

opts_chunk$set(fig.width = 11, fig.height = 10.5, fig.align = 'center') 
opts_chunk$set(tidy=TRUE)
opts_chunk$set(comment = NA)


suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(kableExtra))

strt1<-Sys.time()

#'Loading samples
#'===================

all.samples <- readRDS("03_normalize/allSamplesDataTable.RDS")

print(head(all.samples))

all.samples[,Peptide:= as.character(Peptide),]

setkey(all.samples,Group)

select.samples <- all.samples[J(c("library","AAV","mRNA_53","mRNA_72"))] 

select.samples[,BCcount:=as.integer(mclapply(BC, function(x) length(table(strsplit(paste(t(x), collapse=","), ","))), mc.cores = detectCores()))]
# select.samples[,Score:= BCcount-1,]
select.samples[,Score:= BCcount,]


print("BCcount-1")
print(head(select.samples))

# Check for negative or zero scores and remove such rows
select.samples <- select.samples[Score > 0, ]

select.samples.trsp <- unique(select.samples, by=c("BC","LUTnrs"))

# Check for NA values and remove rows with NA in the Peptide column
select.samples.trsp <- select.samples.trsp[!is.na(Peptide) & Peptide != "", ]

fasta.names <- paste(1:nrow(select.samples.trsp),select.samples.trsp$Score,select.samples.trsp$Group, sep = "|")
write.fasta(as.list(select.samples.trsp$Peptide), fasta.names, "06_motif/trspSamplesPeptides.fasta", open = "w", nbchar = 60, as.string = TRUE)

# # Print the content of the FASTA file for verification
# cat("Content of 06_motif/trspSamplesPeptides.fasta:\n")
# cat(readLines("06_motif/trspSamplesPeptides.fasta"), sep = "\n")

#Generate Scoring table for Weblogo Weighting
select.samples.pepMerge <- select.samples.trsp[, sum(Score), by = c("Peptide")]
setnames(select.samples.pepMerge, "V1", "Score")

print("select.samples.pepMerge")
print(head(select.samples.pepMerge))

#'Executing Hammock Clustering
#'===================

# The Conda environment is named 'HHMER'
conda_env_path <- "~/.conda/envs/HHMER/bin"  # Path to the Conda environment's bin directory

# Set environment variables
Sys.setenv(
  "PATH" = paste(conda_env_path, Sys.getenv("PATH"), sep=":"),
  "HHLIB" = "/home/hooimin/Hammock_v_1.2.0/Hammock_v_1.2.0/hhsuite-2.0.16/lib/hh/"
)

# Example usage of the system command with Hammock
unlink("06_motif/Hammock", recursive = TRUE, force = FALSE)
sys.out <- system(
  paste("java -jar /home/hooimin/Hammock_v_1.2.0/Hammock_v_1.2.0/dist/Hammock.jar full -i 06_motif/trspSamplesPeptides.fasta -d 06_motif/Hammock --max_shift 7 -c 65 -t ", detectCores(), sep = ""),
  intern = TRUE, ignore.stdout = TRUE
)

hammock.log <- data.table(readLines("06_motif/Hammock/run.log"))

print("hammock.log")
print(head(hammock.log))

colnames(hammock.log) <- c("Hammock log file")
knitr::kable(hammock.log, longtable = T)


#'Generation of Weblogo visualization
#'===================
ham.clusters <- data.table(read.table("06_motif/Hammock/final_clusters.tsv", header = TRUE, skip = 0, sep="\t",
                                      stringsAsFactors = FALSE, fill=TRUE))

print("ham.clusters")
print(head(ham.clusters))

id.order <- as.list(ham.clusters$cluster_id)
ham.clusters.all <- data.table(read.table("06_motif/Hammock/final_clusters_sequences.tsv", header = TRUE, skip = 0, sep="\t",
                                          stringsAsFactors = FALSE, fill=TRUE))

print("ham.clusters.all")
print(head(ham.clusters.all))

ham.clusters.all[,alignment := gsub('\\-', '\\_', alignment)]

print("ham.clusters.all after gsub")
print(head(ham.clusters.all))

setkey(select.samples,Peptide)
setkey(select.samples.trsp,Peptide)

print("ham.clusters.all")
print(head(ham.clusters.all))

# Ensure the output directories exist
unlink("06_motif/WEBlogos", recursive = TRUE, force = FALSE)
dir.create(file.path("06_motif/", "WEBlogos"), showWarnings = FALSE)
dir.create(file.path("06_motif/Hammock/", "alignments_final_Scored"), showWarnings = FALSE)

setkey(ham.clusters.all,cluster_id)
setkey(ham.clusters,cluster_id)
setkey(select.samples.pepMerge,Peptide)


opts_chunk$set(out.width='100%', fig.align = 'center')
generateWeblogo <- function(in.name) {
  # Set WebLogo path
  weblogo_env_path <- "~/.conda/envs/weblogo/bin/" 
  # Set Ghostscript path
  ghostscript_path <- "~/.conda/envs/weblogo/bin/" 

  # Update PATH environment variable to include WebLogo Conda environment
  Sys.setenv(
    "PATH" = paste(weblogo_env_path, ghostscript_path, Sys.getenv("PATH"), sep=":")
  )

  # Verify PATH
  print(Sys.getenv("PATH"))
  Sys.setenv(GHOSTSCRIPT = file.path(ghostscript_path, "gs"))

  #in.name <- ham.clusters$cluster_id[2]
  # in.name <- 6777
  
  this.fa <- read.fasta(file = paste("06_motif/Hammock/alignments_final/", in.name, ".aln", sep=""))
  print("this.fa")
  print(head(this.fa))

  allSeqs <- unlist(getSequence(this.fa, as.string = TRUE))

  allSeqs <- data.table(unlist(lapply(allSeqs, function(x) gsub("([-])","",toupper(x)))))

  print("allSeqs")
  print(head(allSeqs))

  allSeqs.out <- select.samples.pepMerge[J(allSeqs)]

  print("allSeqs.out")
  print(head(allSeqs.out))

  allSeqs.out$Annot <- data.table(getName(this.fa))

  print("allSeqs.out1")
  print(head(allSeqs.out))

  allSeqs.out[,Annot:= paste(Annot,"_",Score,sep="")]

  print("allSeqs.out2")
  print(head(allSeqs.out))

  allSeqs.out$Alignment <- data.table(toupper(unlist(getSequence(this.fa, as.string = TRUE))))

  print("allSeqs.out3")
  print(head(allSeqs.out))
  
  allSeqs.out <- allSeqs.out[rep(1:.N,Score)][,Indx:=1:.N,by=Peptide]

  print("allSeqs.out4")
  print(head(allSeqs.out))

  allSeqs.out[,Annot:= paste(Annot,"_",Indx,sep="")]

  print("allSeqs.out5")
  print(head(allSeqs.out))
  
  write.fasta(as.list(allSeqs.out$Alignment), allSeqs.out$Annot, nbchar = 60, paste("06_motif/Hammock/alignments_final_Scored/", in.name, ".aln", sep=""), open = "w")
  
  this.main <- ham.clusters[J(in.name)]

  print("this.main")
  print(head(this.main))

  main.gene <- select.samples.trsp[J(this.main$main_sequence)]$GeneName[1]

  print("main.gene")
  print(head(main.gene))

  # this.title <- paste("## Peptide",this.main$main_sequence,"from",main.gene,"with cluster number",in.name, sep=" ")
  this.title <- paste(this.main$main_sequence, main.gene, "cluster:", in.name, sep=" ")

  
  print("this.title")
  print(head(this.title))

  # Run WebLogo command to generate EPS
  weblogo_cmd_eps <- paste(
    file.path(weblogo_env_path, "weblogo"),
    " --format eps --sequence-type protein --size large --errorbars NO --resolution 300 --composition equiprobable --color-scheme chemistry --title '", 
    this.title, 
    "' < 06_motif/Hammock/alignments_final_Scored/", 
    in.name, 
    ".aln > 06_motif/WEBlogos/", 
    in.name, 
    ".eps", 
    sep = ""
  )
  system(weblogo_cmd_eps, intern = TRUE, ignore.stdout = FALSE)

  # Run Ghostscript command to convert EPS to PDF
  gs_cmd <- paste(
    file.path(ghostscript_path, "gs"),
    " -sDEVICE=pdfwrite -dPDFSETTINGS=/printer -sstdout=%stderr -dColorConversionStrategy=/LeaveColorUnchanged -sOutputFile=06_motif/WEBlogos/", 
    in.name, 
    ".pdf -dDEVICEWIDTHPOINTS=264 -dDEVICEHEIGHTPOINTS=119 -dSAFER -dNOPAUSE -r300 06_motif/WEBlogos/", 
    in.name, 
    ".eps", 
    sep = ""
  )
  system(gs_cmd, intern = TRUE, ignore.stdout = FALSE)

  pdf_file <- paste("06_motif/WEBlogos/", in.name, ".pdf", sep="")
  if (!file.exists(pdf_file)) {
    stop(paste("WebLogo command failed. PDF file not created:", pdf_file))
  }

  cat('\n')
  cat(this.title, "\n")
  cat('\n')
  cat('\n')
  cat(paste0("![Peptide: ", this.main$main_sequence, " from ", main.gene, " with cluster number ", in.name, "](06_motif/WEBlogos/", in.name, ".pdf)"))
  cat('\n')
  
  # Check if this.main has at least 7 columns
  if (ncol(this.main) < 7) {
    stop("this.main does not have enough columns.")
  }
  
  out.table <- knitr::kable(this.main[, c(1:7)], format = "latex")
  print(column_spec(out.table, 1:7, monospace = TRUE) %>% kable_styling(latex_options = c("striped", "scale_down", "repeat_header")))
  cat('\n')
  this.cluster <- ham.clusters.all[J(in.name)]
  out.table <- knitr::kable(this.cluster[, c(1:7)], format = "latex")
  print(column_spec(out.table, 1:7, monospace = TRUE) %>% kable_styling(latex_options = c("striped", "scale_down", "repeat_header")))
  
  this.found <- select.samples[J(this.cluster$sequence)]
  setnames(this.cluster, "sequence", "Peptide")
  this.found <- merge(this.found, this.cluster[, 2:3], by = "Peptide", all = FALSE)
  cat('\n')
  cat('\n')
  output.order <- c("alignment", "LUTnrs", "GeneName", "start", "Group", "Score")
  if (nrow(this.found) >= 48) {
    this.found.p1 <- this.found[1:47,]
    out.table <- knitr::kable(this.found.p1[, ..output.order], format = "latex")
    print(column_spec(out.table, 1, monospace = TRUE) %>% kable_styling(latex_options = c("striped", "repeat_header")))
    cat('\n')
    cat("\n\n\\pagebreak\n")
    cat("\n\n\\clearpage\n")
    this.found <- this.found[48:nrow(this.found),]
    if (nrow(this.found) >= 48) {
      this.found.p2 <- this.found[1:47,]
      out.table <- knitr::kable(this.found.p2[, ..output.order], format = "latex")
      print(column_spec(out.table, 1, monospace = TRUE) %>% kable_styling(latex_options = c("striped", "repeat_header")))
      cat('\n')
      cat("\n\n\\pagebreak\n")
      cat("\n\n\\clearpage\n")
      this.found <- this.found[48:nrow(this.found),]
    }
    if (nrow(this.found) >= 48) {
      this.found.p2 <- this.found[1:47,]
      out.table <- knitr::kable(this.found.p2[, ..output.order], format = "latex")
      print(column_spec(out.table, 1, monospace = TRUE) %>% kable_styling(latex_options = c("striped", "repeat_header")))
      cat('\n')
      cat("\n\n\\pagebreak\n")
      cat("\n\n\\clearpage\n")
      this.found <- this.found[48:nrow(this.found),]
    }
    out.table <- knitr::kable(this.found[, ..output.order], format = "latex")
    print(column_spec(out.table, 1, monospace = TRUE) %>% kable_styling(latex_options = c("striped", "repeat_header")))
    cat('\n')
    cat("\n\n\\clearpage\n")
    cat('\n')
  } else {
    out.table <- knitr::kable(this.found[, ..output.order], format = "latex")
    print(column_spec(out.table, 1, monospace = TRUE) %>% kable_styling(latex_options = c("striped", "scale_down", "repeat_header")))
    cat('\n')
    cat("\n\n\\pagebreak\n")
    cat("\n\n\\clearpage\n")
  }
}

#+ results = 'asis'
invisible(lapply(id.order, generateWeblogo))
#+ results = 'markup'

# setkey(select.samples.trsp, Peptide)
# select.samples.trsp.select <- select.samples.trsp[J(c("PPDELNLTTASLPL"))]

print("Total analysis time:")
print(Sys.time() - strt1)

devtools::session_info()