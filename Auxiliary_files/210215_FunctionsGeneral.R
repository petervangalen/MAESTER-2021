# Peter van Galen, 210215
# General functions for analyses in the MAESTER project

# General
message("cutf()")
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Function that computes all heteroplasmic variants from MAEGATK output (from Caleb Lareau). Rows represents a position along the mitochondrial genome and the three possible disagreements with the reference (except 3107 has four possible disagreements because the reference is N)
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
