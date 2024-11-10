# Clear the workspace
rm(list = ls())

setwd("/home/hooimin/lu2024-17-19/Data_RAAV/p007_data")

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

strt1 <- Sys.time()

#'Loading samples
#'===================
# Define the base directory
base_dir <- "/home/hooimin/lu2024-17-19/Data_RAAV/p007_data"

# Define the path to the 03_normalize subdirectory
normalize_dir <- file.path(base_dir, "03_normalize")
  
# List all subdirectories within the 03_normalize directory
subdirectories <- list.dirs(path = normalize_dir, full.names = TRUE, recursive = TRUE)
  
# Loop through each subdirectory
for (subdir in subdirectories) {
    
    # Extract the base name of the subdirectory
    base_name <- basename(subdir)
    
    # Check if the base_name matches the specified subdirectories
    if (base_name %in% c("38_S13", "42_S15", "44_S17", "49_S19", "51_S21", "43_S16", "46_S18", "50_S20", "52_S22")) {
        # Create a new directory for the current base_name
        output_dir <- file.path("06_motif", base_name)
        dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
        # List RDS files matching the pattern in the current subdirectory
        rds_files <- list.files(path = subdir, pattern = "_normalized_p007_53.*\\.RDS$", full.names = TRUE)
    
        # Debug: Print the subdirectory and RDS files found
        cat("Processing subdirectory:", subdir, "\n")
        cat("RDS files found:", rds_files, "\n")

        # Process each RDS file
        for (rds_file in rds_files) {
            # Read the RDS file
            all.samples <- readRDS(rds_file)
      
            all.samples[, Peptide := as.character(Peptide)]
            setkey(all.samples, Group)
            select.samples <- all.samples[J(c("library", "mRNA_All", "mRNA_53"))] 
            select.samples[, BCcount := as.integer(mclapply(BC, function(x) length(table(strsplit(paste(t(x), collapse = ","), ","))), mc.cores = detectCores()))]
            select.samples[, Score := BCcount]
      
            print("BCcount-1")
            print(head(select.samples))
      
            # Check for negative or zero scores and remove such rows
            select.samples <- select.samples[Score > 0, ]
            select.samples.trsp <- unique(select.samples, by = c("BC", "LUTnrs"))
      
            # Check for NA values and remove rows with NA in the Peptide column
            select.samples.trsp <- select.samples.trsp[!is.na(Peptide) & Peptide != "", ]
      
            fasta.names <- paste(1:nrow(select.samples.trsp), select.samples.trsp$Score, select.samples.trsp$Group, sep = "|")
            fasta_file_path <- file.path(output_dir, paste0(base_name, "_trspSamplesPeptides.fasta"))
            write.fasta(as.list(select.samples.trsp$Peptide), fasta.names, fasta_file_path, open = "w", nbchar = 60, as.string = TRUE)
      
            # Debug: Print the FASTA file path
            cat("FASTA file created:", fasta_file_path, "\n")
      
            # Generate Scoring table for Weblogo Weighting
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
                "PATH" = paste(conda_env_path, Sys.getenv("PATH"), sep = ":"),
                "HHLIB" = "/home/hooimin/Hammock_v_1.2.0/Hammock_v_1.2.0/hhsuite-2.0.16/lib/hh/"
            )
      
            # Example usage of the system command with Hammock
            hammock_output_dir <- file.path(output_dir, "Hammock")
            unlink(hammock_output_dir, recursive = TRUE, force = FALSE)
            hammock_cmd <- paste(
                "java -jar /home/hooimin/Hammock_v_1.2.0/Hammock_v_1.2.0/dist/Hammock.jar full -i", fasta_file_path, "-d", hammock_output_dir, "--max_shift 5 -c 65 -t", detectCores(), "--greedy_threshold 15",
                sep = " "
            )

            # Run the Hammock command and handle errors
            tryCatch({
                sys.out <- system(hammock_cmd, intern = TRUE, ignore.stdout = TRUE)
                
                hammock_log_path <- file.path(hammock_output_dir, "run.log")
                hammock.log <- data.table(readLines(hammock_log_path))
                
                # Check for StringIndexOutOfBoundsException in the log
                if (any(grepl("String index out of range: 15", hammock.log))) {
                    message("String index out of range error detected. Skipping file: ", fasta_file_path)
                    next
                }
                
                print("hammock.log")
                print(head(hammock.log))
                
                colnames(hammock.log) <- c("Hammock log file")
                knitr::kable(hammock.log, longtable = T)
                
                #'Generation of Weblogo visualization
                #'===================
                ham_clusters_path <- file.path(hammock_output_dir, "final_clusters.tsv")
                ham.clusters <- data.table(read.table(ham_clusters_path, header = TRUE, skip = 0, sep = "\t", stringsAsFactors = FALSE, fill = TRUE))
                
                print("ham.clusters")
                print(head(ham.clusters))
                
                id.order <- as.list(ham.clusters$cluster_id)
                ham_clusters_all_path <- file.path(hammock_output_dir, "final_clusters_sequences.tsv")
                ham.clusters.all <- data.table(read.table(ham_clusters_all_path, header = TRUE, skip = 0, sep = "\t", stringsAsFactors = FALSE, fill = TRUE))
                
                print("ham.clusters.all")
                print(head(ham.clusters.all))
                
                ham.clusters.all[, alignment := gsub('\\-', '\\_', alignment)]
                
                print("ham.clusters.all after gsub")
                print(head(ham.clusters.all))
                
                setkey(select.samples, Peptide)
                setkey(select.samples.trsp, Peptide)
                
                print("ham.clusters.all")
                print(head(ham.clusters.all))
                
                # Ensure the output directories exist
                weblogos_dir <- file.path(output_dir, "WEBlogos")
                alignments_final_scored_dir <- file.path(hammock_output_dir, "alignments_final_Scored")
                unlink(weblogos_dir, recursive = TRUE, force = FALSE)
                dir.create(weblogos_dir, showWarnings = FALSE)
                dir.create(alignments_final_scored_dir, showWarnings = FALSE)
                
                setkey(ham.clusters.all, cluster_id)
                setkey(ham.clusters, cluster_id)
                setkey(select.samples.pepMerge, Peptide)
                
                opts_chunk$set(out.width = '100%', fig.align = 'center')
                generateWeblogo <- function(in.name) {
                    # Set WebLogo path
                    weblogo_env_path <- "~/.conda/envs/weblogo/bin/" 
                    # Set Ghostscript path
                    ghostscript_path <- "~/.conda/envs/weblogo/bin/" 
                    
                    # Update PATH environment variable to include WebLogo Conda environment
                    Sys.setenv(
                        "PATH" = paste(weblogo_env_path, ghostscript_path, Sys.getenv("PATH"), sep = ":")
                    )
                    
                    # # Verify PATH
                    # print(Sys.getenv("PATH"))
                    # Sys.setenv(GHOSTSCRIPT = file.path(ghostscript_path, "gs"))
                    
                    this_fa_path <- file.path(hammock_output_dir, "alignments_final", paste0(in.name, ".aln"))
                    this.fa <- read.fasta(this_fa_path)
                    print("this.fa")
                    print(head(this.fa))
                    
                    allSeqs <- unlist(getSequence(this.fa, as.string = TRUE))
                    
                    allSeqs <- data.table(unlist(lapply(allSeqs, function(x) gsub("([-])", "", toupper(x)))))
                    
                    allSeqs.out <- select.samples.pepMerge[J(allSeqs)]
                    
                    allSeqs.out$Annot <- data.table(getName(this.fa))
                    
                    allSeqs.out[, Annot := paste(Annot, "_", Score, sep = "")]
                    
                    allSeqs.out$Alignment <- data.table(toupper(unlist(getSequence(this.fa, as.string = TRUE))))
                    
                    allSeqs.out <- allSeqs.out[rep(1:.N, Score)][, Indx := 1:.N, by = Peptide]
                    
                    allSeqs.out[, Annot := paste(Annot, "_", Indx, sep = "")]
                    
                    write_fasta_path <- file.path(alignments_final_scored_dir, paste0(in.name, ".aln"))
                    write.fasta(as.list(allSeqs.out$Alignment), allSeqs.out$Annot, nbchar = 60, write_fasta_path, open = "w")
                    
                    this.main <- ham.clusters[J(in.name)]
                    
                    # Check the columns of this.main
                    cat("Columns in this.main:", colnames(this.main), "\n")
                    cat("Number of columns in this.main:", ncol(this.main), "\n")
                    
                    if (ncol(this.main) < 6) {
                      cat("Warning: this.main in subdirectory", base_name, "does not have enough columns.\n")
                      return(NULL)
                    }
                    print("this.main")
                    print(head(this.main))
                    
                    main.gene <- select.samples.trsp[J(this.main$main_sequence)]$GeneName[1]
                    
                    print("main.gene")
                    print(head(main.gene))
                    
                    # this.title <- paste("## Peptide", this.main$main_sequence, "from", main.gene, "with cluster number", in.name, sep = " ")
                    this.title <- paste(this.main$main_sequence, main.gene, "cluster:", in.name, sep = " ")
                    
                    print("this.title")
                    print(head(this.title))
                    
                    # Run WebLogo command to generate EPS
                    weblogo_cmd_eps <- paste(
                      file.path(weblogo_env_path, "weblogo"),
                      " --format eps --sequence-type protein --size large --errorbars NO --resolution 300 --composition equiprobable --color-scheme chemistry --title '", 
                      this.title, 
                      "' < ", write_fasta_path, " > ", file.path(weblogos_dir, paste0(in.name, ".eps")),
                      sep = ""
                    )
                    system(weblogo_cmd_eps, intern = TRUE, ignore.stdout = FALSE)
                    
                    # Run Ghostscript command to convert EPS to PDF
                    gs_cmd <- paste(
                      file.path(ghostscript_path, "gs"),
                      " -sDEVICE=pdfwrite -dPDFSETTINGS=/printer -sstdout=%stderr -dColorConversionStrategy=/LeaveColorUnchanged -sOutputFile=", file.path(weblogos_dir, paste0(in.name, ".pdf")),
                      " -dDEVICEWIDTHPOINTS=264 -dDEVICEHEIGHTPOINTS=119 -dSAFER -dNOPAUSE -r300 ", file.path(weblogos_dir, paste0(in.name, ".eps")),
                      sep = ""
                    )
                    system(gs_cmd, intern = TRUE, ignore.stdout = FALSE)
                    
                    pdf_file <- file.path(weblogos_dir, paste0(in.name, ".pdf"))
                    if (!file.exists(pdf_file)) {
                      stop(paste("WebLogo command failed. PDF file not created:", pdf_file))
                    }
                    
                    cat('\n')
                    cat(this.title, "\n")
                    cat('\n')
                    cat('\n')
                    cat(paste0("![Peptide: ", this.main$main_sequence, " from ", main.gene, " with cluster number ", in.name, "](06_motif/WEBlogos/", in.name, ".pdf)"))
                    cat('\n')
                    
                    out.table <- knitr::kable(this.main[, c(1:6)], format = "latex")
                    print(column_spec(out.table, 1:7, monospace = TRUE) %>% kable_styling(latex_options = c("striped", "scale_down", "repeat_header")))
                    cat('\n')
                    this.cluster <- ham.clusters.all[J(in.name)]
                    out.table <- knitr::kable(this.cluster[, c(1:6)], format = "latex")
                    print(column_spec(out.table, 1:6, monospace = TRUE) %>% kable_styling(latex_options = c("striped", "scale_down", "repeat_header")))
                    
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
      
                print("Total analysis time:")
                print(Sys.time() - strt1)
      
                devtools::session_info()
            }, error = function(e) {
                message("Error processing sample ", base_name, ": ", e$message)
                # Skip to the next iteration if an error occurs
            })
        }
    }
}

# Run the Java program to check the length of the string
java_output <- system("java HHsuiteRunner", intern = TRUE)
cat(java_output, "\n")