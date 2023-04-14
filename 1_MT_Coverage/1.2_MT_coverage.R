# Peter van Galen, 210924
# Plot coverage of chrM from glioblastoma maegatk data objects

# Prerequisites -----------------------------------------------------------------------------------

options(scipen = 999)

library(tidyverse)
library(Matrix)
library(ggforce)
library(SummarizedExperiment)

rm(list=ls())
setwd("~/DropboxPartners/Projects/Maester/AnalysisPeter/1_MT_Coverage")

# Functions (available at https://github.com/vangalenlab/MAESTER-2021)
source("../210215_FunctionsGeneral.R")


# Load locally saved maegatk data (choose one) ----------------------------------------------------
# This data is available at https://vangalenlab.bwh.harvard.edu/maester-2021/
    experiment.name <- "SW_MGH252_A"
    maegatk.full <- readRDS(file = "SW_MGH252_A_mr3_maegatk.rds")
    metadata.df <- read.table("../8_Glioblastoma/Chadi email 210721/MGH252_A_C_PBMC_MetaData.txt")
    cellMerge4 <- rownames(metadata.df)
    
    experiment.name <- "SW_MGH252_C"
    maegatk.full <- readRDS(file = "SW_MGH252_C_mr3_maegatk.rds")
    metadata.df <- read.table("../8_Glioblastoma/Chadi email 210721/MGH252_A_C_PBMC_MetaData.txt")
    cellMerge4 <- rownames(metadata.df)
    
    experiment.name <- "SW_MGH252_PBMC"
    maegatk.full <- readRDS(file = "SW_MGH252_PBMC_mr3_maegatk.rds")
    metadata.df <- read.table("../8_Glioblastoma/Chadi email 210721/MGH252_A_C_PBMC_MetaData.txt")
    cellMerge4 <- rownames(metadata.df)
    

# Use common cell barcodes for RNAseq and maegatk -------------------------------------------------
# This is not necessary if you subsetted the bam file for high-quality cell barcodes when runnig maegatk
# Only keep cells with a cellMerge id that occurs once, intersect, plot
cellMerge3 <- tibble(cell = cellMerge4) %>% group_by(cell) %>% filter(n()==1) %>% .$cell %>% cutf(d = "\\.", f = 2)
cellMerge2 <- str_c(cellMerge3, "-1")
cellMerge1 <- intersect(cellMerge2, colnames(maegatk.full))
barplot(c(length(cellMerge3), length(cellMerge2), length(cellMerge1)), ylim = c(0, length(cellMerge3)), ylab = "cell number")
axis(side = 1, at = c(0.7, 1.9, 3.1), labels = c("RNAseq_all", "occur_once", "common"))
# Subset maegatk data for high quality cells
maegatk <- maegatk.full[,cellMerge1]
ncol(maegatk.full); ncol(maegatk)


pdf(file = paste0("1.2_", experiment.name, ".pdf"))


# Plot coverage per position ----------------------------------------------------------------------
# Set y axis parameters
ymax <- 200

# Gene locations
GenePos.tib <- tibble(Names = c("MT.ATP6", "MT.ATP8", "MT.CO1", "MT.CO2", "MT.CO3", "MT.CYB", "MT.ND1", "MT.ND2", "MT.ND3",
                                "MT.ND4", "MT.ND4L", "MT.ND5", "MT.ND6", "MT.RNR1", "MT.RNR2"),
                      start = c(8527, 8366, 5904, 7586, 9207, 14747, 3307, 4470, 10059, 10760, 10470, 12337, 14149, 648, 1671), 
                      end = c(9207, 8572, 7445, 8269, 9990, 15887, 4262, 5511, 10404, 12137, 10766, 14148, 14673, 1601, 3229))
GenePos.tib <- GenePos.tib %>% arrange(start) %>%
    mutate(mid = round((end-start)/2+start,0), ycoord = rep(c(ymax*1.2,ymax*1.1), length.out = 15))

# Plot coverage along chrM
base.tib <- tibble(base = 1:16569, depth = rowMeans(assays(maegatk)[["coverage"]]))

print(
    base.tib %>% ggplot() +
        geom_bar(aes(x = base, y = ifelse(depth > 1, yes = depth, no = NA)), stat = "identity", fill = "#64b53b", width = 1) + 
        coord_cartesian(ylim = c(1, ymax), xlim = c(700, 15900)) +
        scale_y_continuous(trans = "log10") +
        geom_segment(data = GenePos.tib, aes(x = start, y = ycoord, xend = end, yend = ycoord)) +
        geom_text(data = GenePos.tib, aes(x = mid, y = ycoord-ymax*0.2, label = cutf(Names, d = "\\.", f = 2)), size = 3) +
        ylab("Mean coverage per cell") + xlab("Position along chrM") +
        theme_classic() +
        theme(aspect.ratio = 0.5)
)


# Plot mean coverage for top 500 cells ------------------------------------------------------------
cells.tib <- tibble(cell = colnames(maegatk),
                    depth = maegatk$depth)
topcells.tib <- cells.tib %>% slice_max(order_by = depth, n = 500)

print(
    ggplot(topcells.tib, aes(x = 1, y = depth)) +
        geom_violin() +
        geom_sina(size = 0.3) +
        coord_cartesian(ylim = c(0.1, 800)) +
        scale_y_continuous(trans = "log10") +
        ylab("Mean coverage per cell") + xlab("") +
        annotate("text", x = 1, y = max(topcells.tib$depth)*1.5,
                 label = round(mean(topcells.tib$depth), 2)) +
        theme_classic() +
        theme(aspect.ratio = 2, plot.title = element_text(hjust = 0.5)) +
        ggtitle("Mean coverage of top 500 cells")
)


# Plot mean depth for top 5000 bases --------------------------------------------------------------
top.tib <- base.tib %>% arrange(desc(depth)) %>% mutate(key = row_number(), .before = 1) %>% filter(key %in% 1:5000)
print(
    ggplot(top.tib) +
        geom_bar(aes(x = key, y = depth), stat = "identity", fill = "#64b53b", width = 1) +
        coord_cartesian(ylim = c(1, ymax)) +
        scale_y_continuous(trans = "log10") +
        geom_label(data = data.frame(), aes(x = 2500, y = mean(top.tib$depth), label = round(mean(top.tib$depth), 2)),
                   fill = "#64b53b") +
        ylab("Mean coverage per cell") + xlab("Rank sorted position") +
        ggtitle("Mean coverage of top 5000 bases") +
        theme_classic() +
        theme(aspect.ratio = 2)
)

dev.off()



