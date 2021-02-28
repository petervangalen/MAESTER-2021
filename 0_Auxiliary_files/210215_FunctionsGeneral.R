# Peter van Galen, 210215
# General functions for analyses in the Maester project

#### General ####
message("cutf()")
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

#### Maegtk function that computes all heteroplasmic variants from output (from Caleb Lareau) ####
message("computeAFMutMatrix()")
computeAFMutMatrix <- function(SE){
  cov <- assays(SE)[["coverage"]]+ 0.000001
  ref_allele <- as.character(rowRanges(SE)$refAllele)
  
  getMutMatrix <- function(letter){
    mat <- (assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / cov
    rownames(mat) <- paste0(as.character(1:dim(mat)[1]), "_", toupper(ref_allele), ">", letter)
    return(mat[toupper(ref_allele) != letter,])
  }
  
  rbind(getMutMatrix("A"), getMutMatrix("C"), getMutMatrix("G"), getMutMatrix("T"))
}
