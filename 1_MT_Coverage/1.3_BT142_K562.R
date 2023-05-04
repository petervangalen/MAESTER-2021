# Peter van Galen, 210715
# Plot coverage along chrM split by cell lines for reviewer 2 comment 2

# Prerequisites -----------------------------------------------------------------------------------
options(scipen = 999)

library(tidyverse)
library(ggrastr)
library(ggrepel)
library(Matrix)
library(SummarizedExperiment)
library(Seurat)

setwd("~/DropboxMGB/Projects/Maester/AnalysisPeter/1_MT_Coverage/")

rm(list=ls())

# Functions (available at https://github.com/vangalenlab/MAESTER-2021)
source("../Auxiliary_files/210215_FunctionsGeneral.R")


# Load data ---------------------------------------------------------------------------------------

# Locally saved SW maegatk data (available at https://vangalenlab.bwh.harvard.edu/maester-2021/)
tech <- "Seq-Well"
all.maegatk <- readRDS(file = "../1_MT_Coverage/SW_CellLineMix_All_mr3_maegatk.rds")
cells.tib <- readRDS("../3_Cell_line_mixes_variants/SW_CellLineMix_Cells.rds")
bt142.maegatk <- all.maegatk[,filter(cells.tib, CellType_MT == "BT142")$cell]
k562.maegatk <- all.maegatk[,filter(cells.tib, CellType_MT == "K562")$cell]

# Or 10x maegatk data
tech <- "TenX"
all.maegatk <- readRDS(file = "../1_MT_Coverage/TenX_CellLineMix_All_mr3_maegatk.rds")
#seu <- readRDS("../2_Cell_line_mixes_decontX/TenX_CellLineMix_Seurat_Keep.rds")
cells.tib <- readRDS("../3_Cell_line_mixes_variants/TenX_CellLineMix_Cells.rds")
bt142.maegatk <- all.maegatk[,filter(cells.tib, CellType_MT == "BT142")$cell]
k562.maegatk <- all.maegatk[,filter(cells.tib, CellType_MT == "K562")$cell]


# Plot coverage per position ----------------------------------------------------------------------
# Set y axis parameters
ymax <- 1000

# Gene locations
GenePos.tib <- tibble(Names = c("MT.ATP6", "MT.ATP8", "MT.CO1", "MT.CO2", "MT.CO3", "MT.CYB", "MT.ND1", "MT.ND2", "MT.ND3",
                                "MT.ND4", "MT.ND4L", "MT.ND5", "MT.ND6", "MT.RNR1", "MT.RNR2"),
                      start = c(8527, 8366, 5904, 7586, 9207, 14747, 3307, 4470, 10059, 10760, 10470, 12337, 14149, 648, 1671), 
                      end = c(9207, 8572, 7445, 8269, 9990, 15887, 4262, 5511, 10404, 12137, 10766, 14148, 14673, 1601, 3229))
GenePos.tib <- GenePos.tib %>% arrange(start) %>%
    mutate(mid = round((end-start)/2+start,0), ycoord = rep(c(ymax*1.2,ymax*1.1), length.out = 15))

# Plot
base.tib <- tibble(base = 1:16569,
                   bt142_depth = rowMeans(assays(bt142.maegatk)[["coverage"]]),
                   k562_depth = rowMeans(assays(k562.maegatk)[["coverage"]]))

pdf(str_c("1.3_", tech, "_BT142-K562.pdf"))
print(
    base.tib %>% ggplot() +
        geom_bar(aes(x = base, y = ifelse(k562_depth > 1, yes = k562_depth, no = NA)), stat = "identity", fill = "#7BF581", width = 1) + 
        geom_bar(aes(x = base, y = ifelse(bt142_depth > 1, yes = bt142_depth, no = NA)), stat = "identity", fill = "#8E87F5", width = 1) +
        coord_cartesian(ylim = c(1, ymax), xlim = c(700, 15900)) +
        scale_y_continuous(trans = "log10") +
        geom_segment(data = GenePos.tib, aes(x = start, y = ycoord, xend = end, yend = ycoord)) +
        geom_text(data = GenePos.tib, aes(x = mid, y = ycoord-ymax*0.2, label = cutf(Names, d = "\\.", f = 2)), size = 3) +
        ylab("Mean coverage per cell") + xlab("Position along chrM") +
        theme_classic() +
        theme(aspect.ratio = 0.5)
)
dev.off()

