# Peter van Galen, 200110
# Goal: append cell barcode and unique molecular identifier, _CB_UMI, from Read 1 to each Read 2 identifier

options(max.print = 500)
options(stringsAsFactors = FALSE)
options(scipen = 999)


### Load libraries
suppressMessages(library(ShortRead))

rm(list=ls())

### Functions
cutf <- function(x, f=1, d="/", ...) sapply(strsplit(x, d), function(i) i[f], ...)


### Arguments to be provided when executing script
Folder <- commandArgs(trailingOnly=TRUE)[1] # One or multiple directories containing fastq files, not searched recursively
SampleName <- commandArgs(trailingOnly=TRUE)[2] # Sample name that will be used for output files
CellBarcodes <- commandArgs(trailingOnly=TRUE)[3] # Allowlist of cell barcodes to filter by. Could be cells from the CellRanger filtered_feature_bc_matrix, or that passed scRNA-seq QC, or all whitelisted 10x cell barcodes.
CBlength <- as.numeric( commandArgs(trailingOnly=TRUE)[4] ) # Length of the cell barcode (16 for 10x 3' v3)
UMIlength <- as.numeric( commandArgs(trailingOnly=TRUE)[5] ) # Length of the UMI (12 for 10x 3' v3)

### Example values
#Folder <- "/broad/vangalenlab/vangalen/Fastq/201112_NovaSeq/Maester"
#SampleName <- "CellLineMix_Maester"
#CellBarcodes <- "/broad/vangalenlab/vangalen/Maester/1_AssembleFastq/200902/CellLineMix_10X_Allowlist.txt"
#CBlength <- 16
#UMIlength <- 12

### Find R1 fastq files (R2 is found by substitution later)
R1.ch <- list.files(Folder, pattern = paste0(SampleName, ".*_R1_.*fastq.gz$"), full.names = T)
message(Sys.time(), "\nLoading ", length(R1.ch)*2, " fastq files:")
message(cat(c(R1.ch, sub("_R1_", "_R2_", R1.ch)), sep = "\n"))
if(length(R1.ch) == 0) stop("Did not find fastq files.")


### Filter for cell barcodes in the allowlist. Remove -1 from the end (if added by CellRanger count).
cells.split <- unlist( strsplit(CellBarcodes, ",") )
cells.df <- do.call(rbind, lapply(cells.split, function(x) read.table(x)))
cells.ch <- cutf(cells.df$V1, d = "-", f = 1)
if(length(cells.ch) == 0) stop("No cells found.")
message("Found ", length(cells.ch), " cells.\n")


### Process fastq files
message("Read R1 and R2 sequences, filter by ", length(cells.ch), " cell barcodes, write assembled fastq...")
report.ls <- list()  # empty list to store number of reads

# For each R1 fastq
for(f1 in R1.ch) {
  #f1 <- R1.ch[1]
  # Identify R2 fastq
  f2 <- sub("_R1_", "_R2_", f1)

  # Load file in 1E7 read increments
  message("file ", match(f1, R1.ch), "/", length(R1.ch), ": ", basename(f1), " ", appendLF = FALSE)
  strm1 <- FastqStreamer(f1, n=1E7)  # 1M reads by default
  strm2 <- FastqStreamer(f2, n=1E7)

  # For every 10 million reads...  
  repeat{
    message("*", appendLF = FALSE)
    fq1 <- yield(strm1)
    fq2 <- yield(strm2)
    if(length(fq1) == 0 | length(fq2) == 0) break

    # Match to expected cell barcodes
    fq1.m <- ifelse(is.element(as.vector(subseq(sread(fq1), 1, CBlength)), cells.ch), yes = T, no = F)
    report.ls[[basename(f1)]] <- c(length(fq1.m), sum(fq1.m))

    # Filter unmatched reads from the ShortRead objects
    fq1.f <- fq1[fq1.m]
    fq2.f <- fq2[fq1.m]

    # Extract cell barcode and umi from Read1
    fq1.f.cell <- as.vector(subseq(sread(fq1.f), 1, CBlength))
    fq1.f.umi <- as.vector(subseq(sread(fq1.f), CBlength+1, CBlength+UMIlength))

    # Add cell barcode and umi to id of Read2
    #fq2.f@id <- BStringSet(paste(as.vector(fq2.f@id), fq1.f.cell, fq1.f.umi))
    fq2.f@id <- BStringSet(paste0(sub(" .:N:0:", "_", as.vector(fq2.f@id)), "_", fq1.f.cell, "_", fq1.f.umi))
    
    # Check if all the ids of Read1 and Read2 files match up
    if(! all(cutf(as.vector(fq1.f@id), d = " ", f = 1) == cutf(as.vector(fq2.f@id), d = "_", f = 1))) stop("Read ID mismatch")

    # Save fastq file
    writeFastq(fq2.f, file = paste0(SampleName, ".fastq.gz"), mode = "a")
  
    # Clean up memory
    invisible(gc())
  }

  close(strm1)
  close(strm2)
  message(" done")
  invisible(gc())
}


### Generate and save report. NOTE: report.ls only contains data from the last iteration of the repeat{} loop above.
report.mat <- do.call(rbind, report.ls)
report.mat <- rbind(report.mat, colSums(report.mat))
rownames(report.mat)[nrow(report.mat)] <- "total"
report.mat <- cbind(report.mat, report.mat[,2] / report.mat[,1])
colnames(report.mat) <- c("all", "filtered", "fraction")
write.table(report.mat, file = paste0(SampleName, ".stats.txt"), sep = "\t", quote = F)

invisible(gc())

message("\nMaintained ", round(report.mat["total", "fraction"]*100, 2), "% of reads.\n")
sessionInfo()
message("\nFinished!")
