#' ---
#' title: "Converting Fragments Final to index file"
#' author: "Hooi Min Tan Grahn"

(library(data.table))


#'Load sequences
#'===================
LUT.dna <- fread("RAAV_2_mRNA_fragmentation.csv")

#'Remove constitutive backbone sequences
#'===================
LUT.dna[, Seq := gsub("aacctccagagaggcaacgct","",Seq)]
LUT.dna[, Seq := gsub("gccagacaagcagctaccgca","",Seq)]
LUT.dna[, Seq := toupper(Seq)]
setkey(LUT.dna, "Seq")
LUT.dna <- unique(LUT.dna)
LUT.dna[, Names := Seq]
LUT.dna[, LUTnr := make.names(seq(.N), unique=TRUE)]


# Save the processed data to a CSV file
fwrite(LUT.dna, file = "LUTdna.csv")

save(LUT.dna,file = "LUTdna.rda")


