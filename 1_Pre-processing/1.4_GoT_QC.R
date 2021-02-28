# Daniel Ssozi, 201116
# Goal: import the summary table from IronThrone-GoT software, perform QC (ReadThreshold) and save cells for follow-up analyses
# Note, this script only works with R version >= 3.6

options(stringsAsFactors = FALSE)

library(data.table)

rm(list=ls())

# Functions
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Adjustable parameters. Examples are commented out.
dir.ch <- commandArgs(trailingOnly=TRUE)[1]
summtable <- commandArgs(trailingOnly=TRUE)[2]
sample.ch <- commandArgs(trailingOnly=TRUE)[3]
gene.name <-  commandArgs(trailingOnly=TRUE)[4]
# Example values:
#dir.ch <- "201109_BPDCN712_BPDCN712R/ZNF750.1562/"
#summtable <- "201109_BPDCN712_BPDCN712R/ZNF750.1562/ZNF750.1562.summTable.txt"
#sample.ch <- "BPDCN712"
#gene.name <- "ZNF750"

ReadThreshold <- 3
dir.create(dir.ch)


#====================================
# Separate data by UMIs (not cells)
#====================================

# Read data
summtable.df <- read.table(summtable, header = T)
#View(summtable.df)

# Remove GoT calls so we can make our own
summtable.red.df <- summtable.df[,1:(ncol(summtable.df)-3)]

# Make a list of cells
summtable.ls <- split(summtable.red.df, seq(nrow(summtable.red.df)))

# Check if there are always the same number of semi-colon separated values for each cell. This should be TRUE!
stopifnot( all(unlist( lapply(summtable.ls, function(x) length( unique( apply(x, 2, function(x) length(unlist(strsplit(x, split = ";")))))))) == 1) )

# Make a new list where UMIs are separated into rows
UMI.rows.ls <- lapply(summtable.ls, function(x) {
  #x <- summtable.ls[[1]]
  umi.number <- length(unlist( strsplit(x$UMI, split = ";") ))
  current.cell <- do.call(rbind, lapply(1:umi.number, function(y) apply(x, 2, function(z) cutf(z, d = ";", f = y)) ))
  data.frame( current.cell )
} )

# Make data frame with one UMI per row (instead of one cell)
summtable.UMI.df <- do.call(rbind, UMI.rows.ls)
# Make some columns (that we use later) numeric
summtable.UMI.df$num.WT.in.dups <- as.numeric(summtable.UMI.df$num.WT.in.dups)
summtable.UMI.df$num.MUT.in.dups <- as.numeric(summtable.UMI.df$num.MUT.in.dups)


#=============================
# Subset UMI Summary Table
#=============================

# Split the summary table into two data frames (wt/mut) with rows indicating cell barode, UMI, and the number of reads
summtable.UMI.wt.df <- summtable.UMI.df[summtable.UMI.df$num.WT.in.dups > 0, c("BC", "UMI", "num.WT.in.dups")]
colnames(summtable.UMI.wt.df) <- c("BC", "UMI", "reads")
summtable.UMI.mut.df <- summtable.UMI.df[summtable.UMI.df$num.MUT.in.dups > 0,c("BC", "UMI", "num.MUT.in.dups")]
colnames(summtable.UMI.mut.df) <- c("BC", "UMI", "reads")

# Similarly, split the summary table into two data frames, but only with high quality calls using the ReadThreshold
summtable.UMI.wt.subset.df <- summtable.UMI.df[summtable.UMI.df$num.WT.in.dups >= ReadThreshold & summtable.UMI.df$num.WT.in.dups >= ReadThreshold*summtable.UMI.df$num.MUT.in.dups, c("BC", "UMI", "num.WT.in.dups")]
colnames(summtable.UMI.wt.subset.df) <- c("BC", "UMI", "reads")
summtable.UMI.mut.subset.df <- summtable.UMI.df[summtable.UMI.df$num.MUT.in.dups >= ReadThreshold & summtable.UMI.df$num.MUT.in.dups >= ReadThreshold*summtable.UMI.df$num.WT.in.dups, c("BC", "UMI", "num.MUT.in.dups")]
colnames(summtable.UMI.mut.subset.df) <- c("BC", "UMI", "reads")


#=============================
# Line plot
#=============================

# Generate line plot of reads for each UMI, sorted from the most to the least reads, lines separated by wt or mut reads

# Determine limits for axes scaling
max.transcripts <- max(nrow(summtable.UMI.wt.df), nrow(summtable.UMI.mut.df))
max.reads <- max(summtable.UMI.wt.df$reads, summtable.UMI.mut.df$reads)

# Plot
pdf(file = paste0(dir.ch, gene.name, ".QC.pdf"), width = 12, height = 12)
par(mar=c(4,4,4,4), pty = "s", par(mfrow=c(2,2)))

plot(NA, xlab = "Transcripts", ylab = "Reads", log = "y", xlim=c(1, max.transcripts), ylim=c(1, max.reads), main = paste0("Reads per transcript in ", sample.ch))
abline(h=ReadThreshold, lty=2)
lines(sort(summtable.UMI.wt.df$reads, decreasing = T), col="black")
lines(sort(summtable.UMI.mut.df$reads, decreasing = T), col="red")
text(x = max.transcripts, y = max.reads*0.6, labels = paste0(gene.name, "_wt"), pos=2)
text(x = max.transcripts, y = max.reads*0.4, labels = paste0(gene.name, "_mut"), pos=2, col="red")


#====================================
# Scatter plot
#====================================

# Plot wild-type vs mutant reads for each transcript to indicate consistency of wt / mut calls between reads
plot(x = summtable.UMI.df$num.WT.in.dups+1, y = summtable.UMI.df$num.MUT.in.dups+1,
     pch = 1, 
     log = "xy", xlab = "Wild-type calls + 1", ylab = "Mutant calls + 1", main = "Reads per transcript: Wt vs Mut",
     xlim = c(1, max(100, summtable.UMI.df$num.WT.in.dups+1)), 
     ylim = c(1, max(100, summtable.UMI.df$num.MUT.in.dups+1)))

# Add diagonal lines to indicate the QC parameters
abline(log10(ReadThreshold), 1, lty = 2)
abline(-log10(ReadThreshold), 1, lty = 2)

# Add numbers to indicate how many dots fall in the upper left or lower right parts of the graph
text(x = max(100, summtable.UMI.wt.subset.df$reads+1)*0.9, y = 2, labels = nrow(summtable.UMI.wt.subset.df))
text(x = 2, y = max(100, summtable.UMI.mut.subset.df$reads+1)*0.9, labels = nrow(summtable.UMI.mut.subset.df))


#============================================
# Make a data table with information per cell 
#============================================

# Make a list of data.frames, one element for each cell. 
wt.filter.cells.ls <- split(summtable.UMI.wt.subset.df, summtable.UMI.wt.subset.df$BC)
# Make a data.table where each row represents a cell, and columns indicate the number of wtUMIs (number of rows in each list element) and wtReads (sum of the reads column in each list element)
wt.UMIs.per.cell.dt <- data.table(BC = names(wt.filter.cells.ls),
                                  wtUMIs = 0 + as.numeric( do.call(rbind, lapply(wt.filter.cells.ls, nrow)) ),
                                  wtReads = 0 + unlist(lapply(wt.filter.cells.ls, function(x) sum(x$reads))))
setorder(wt.UMIs.per.cell.dt, -wtUMIs)

### Do the same for mutant
mut.filter.cells.ls <- split(summtable.UMI.mut.subset.df, summtable.UMI.mut.subset.df$BC)
mut.UMIs.per.cell.dt <- data.table(BC = names(mut.filter.cells.ls),
                                   mutUMIs = 0 + as.numeric( do.call(rbind, lapply(mut.filter.cells.ls, nrow))),
                                   mutReads = 0 + unlist(lapply(mut.filter.cells.ls, function(x) sum(x$reads))))
setorder(mut.UMIs.per.cell.dt, -mutUMIs)

# Merge wt and mut, order from most to least cells
Cell_info.dt <- merge(x = wt.UMIs.per.cell.dt, y = mut.UMIs.per.cell.dt, by = "BC", all = T)
Cell_info.dt[is.na(Cell_info.dt)] <- 0
Cell_info.dt <- Cell_info.dt[order(rowSums(Cell_info.dt[,c("wtReads", "mutReads")]), decreasing = T),]

if(nrow(Cell_info.dt) == 0) {
  stop("No wild-type or mutant transcripts passed QC thresholds.")
}


#================================
# Bar plots
#================================

# Make a bar plot indicating the number of reads per cell (stacked for wt/mut)
Read.matrix.mx <- data.matrix(Cell_info.dt[,c("wtReads", "mutReads")])
rownames(Read.matrix.mx) <- Cell_info.dt$BC
barplot(t(Read.matrix.mx),
        las = 2,
        main = "Reads per cell",
        ylab = "total reads",
        cex.names = 1/3,
        col = c("grey", "red"),
        border = NA)

# Make a bar plot indicating the number of UMIs per cell (stacked for wt/mut)
UMI.matrix.mx <- data.matrix(Cell_info.dt[,c("wtUMIs", "mutUMIs")])
rownames(UMI.matrix.mx) <- Cell_info.dt$BC
barplot(t(UMI.matrix.mx),
        las = 2,
        main = "UMI per cell",
        ylab = "total UMIs",
        cex.names = 1/3,
        col = c("grey", "red"),
        border = NA)

dev.off()


#================================
# Save results
#================================


write.table(Cell_info.dt[,c("BC", "wtUMIs", "mutUMIs")], file = gsub("summTable", "FilteredCells", summtable), quote = F, row.names = F, sep = "\t")
sessionInfo()

