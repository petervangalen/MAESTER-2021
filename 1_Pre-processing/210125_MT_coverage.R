# Peter van Galen, 210125
# Plot coverage of chrM from maegatk data object

options(stringsAsFactors = FALSE)
options(scipen = 999)

library(tidyverse)
library(Matrix)
library(ggforce)
library(SummarizedExperiment)
library(ggbeeswarm)
#library(data.table)
#library(Seurat)

rm(list=ls())
setwd("~/DropboxPartners/Projects/Maester/AnalysisPeter/1_MT_Coverage")

# More functions (available on https://github.com/vangalenlab/MAESTER-2021)
source("../210215_FunctionsGeneral.R")

# Choose one
experiment.name <- "SW_CellLineMix"
experiment.name <- "TenX_CellLineMix_Combined"
experiment.name <- "TenX_BPDCN712"

# Load locally saved maegatk data.
# This data is available from https://vangalenlab.bwh.harvard.edu/maester-2021/
# Note that these are lists of two maegatk objects (one for scRNA-seq coverage alone, one for scRNA-seq+MAESTER coverage)
se.ls <- readRDS(file = paste0(experiment.name, "_mr3_maegatk.rds"))

pdf(file = paste0(experiment.name, ".pdf"))

#### Plot mean coverage for top 500 cells
# Use same cells for RNAseq and Maester data
common.cells <- intersect(colnames(se.ls[[1]]), colnames(se.ls[[2]]))
se.ls[[1]] <- se.ls[[1]][,common.cells]
se.ls[[2]] <- se.ls[[2]][,common.cells]

# Mean coverage per cell
cells.tib <- tibble(cell = colnames(se.ls[[1]]),
                    rnaseq_depth = se.ls[[1]]$depth,
                    maester_depth = se.ls[[2]]$depth)
topcells.tib <- cells.tib %>% top_n(500)

print(
pivot_longer(topcells.tib, cols = c("rnaseq_depth", "maester_depth")) %>%
    mutate(name = factor(name, levels = c("rnaseq_depth", "maester_depth"))) %>%
    ggplot(aes(x = name, y = value)) +
    geom_violin() +
    geom_sina(size = 0.3) +
    coord_cartesian(ylim = c(0.1, 800)) +
    scale_y_continuous(trans = "log10") +
    ylab("Mean coverage per cell") + xlab("") +
    annotate("text", x = c(1,2), y = c(max(topcells.tib$rnaseq_depth),max(topcells.tib$maester_depth))*1.5,
             label = round(c(mean(topcells.tib$rnaseq_depth), mean(topcells.tib$maester_depth)), 2)) +
    theme_classic() +
    theme(aspect.ratio = 2, plot.title = element_text(hjust = 0.5)) +
    ggtitle("Mean coverage of top 500 cells")
)

### Plot coverage per position
# Gene locations
GenePos.tib <- tibble(Names = c("MT.ATP6", "MT.ATP8", "MT.COX1", "MT.COX2", "MT.COX3", "MT.CYTB", "MT.ND1", "MT.ND2", "MT.ND3",
                                "MT.ND4", "MT.ND4L", "MT.ND5", "MT.ND6", "MT.RNR1", "MT.RNR2"),
                      start = c(8527, 8366, 5904, 7586, 9207, 14747, 3307, 4470, 10059, 10760, 10470, 12337, 14149, 648, 1671), 
                      end = c(9207, 8572, 7445, 8269, 9990, 15887, 4262, 5511, 10404, 12137, 10766, 14148, 14673, 1601, 3229))
GenePos.tib <- GenePos.tib %>% arrange(start) %>%
    mutate(mid = round((end-start)/2+start,0), ycoord = rep(c(250,230), length.out = 15))

# Plot
base.tib <- tibble(base = 1:16569,
                   rnaseq_depth = rowMeans(assays(se.ls[[1]])[["coverage"]]),
                   maester_depth = rowMeans(assays(se.ls[[2]])[["coverage"]]))
print(
base.tib %>% ggplot() +
    geom_bar(aes(x = base, y = ifelse(maester_depth > 1, yes = maester_depth, no = NA)), stat = "identity", fill = "#64b53b", width = 1) + 
    geom_bar(aes(x = base, y = ifelse(rnaseq_depth > 1, yes = rnaseq_depth, no = NA)), stat = "identity", fill = "#fdcb25", width = 1) +
    coord_cartesian(ylim = c(1, 1000)) + # previously 200
    scale_y_continuous(trans = "log10") + #, breaks = c(1, 2, 5, 10, 20, 50, 100, 200)) +
    geom_segment(data = GenePos.tib, aes(x = start, y = ycoord, xend = end, yend = ycoord)) +
    geom_text(data = GenePos.tib, aes(x = mid, y = ycoord-50, label = cutf(Names, d = "\\.", f = 2)), size = 3) +
    ylab("Mean coverage per cell") + xlab("Position along chrM") +
    theme_classic() +
    theme(aspect.ratio = 0.5)
)

#### Plot mean depth for top 5000 bases
top.tib <- base.tib %>% arrange(desc(maester_depth)) %>% mutate(key = row_number(), .before = 1) %>% filter(key %in% 1:5000)
print(
ggplot(top.tib) +
    geom_bar(aes(x = key, y = maester_depth), stat = "identity", fill = "#64b53b", width = 1) +
    geom_bar(aes(x = key, y = ifelse(rnaseq_depth > 1, yes = rnaseq_depth, no = NA)), stat = "identity", fill = "#fdcb25", width = 1) +
    coord_cartesian(ylim = c(1, 1000)) + # previously 200
    scale_y_continuous(trans = "log10") + #, breaks = c(1, 2, 5, 10, 20, 50, 100, 200)) +
    geom_label(data = data.frame(), aes(x = 2500, y = mean(top.tib$maester_depth), label = round(mean(top.tib$maester_depth), 2)),
               fill = "#64b53b") +
    geom_label(data = data.frame(), aes(x = 2500, y = 2, label = round(mean(top.tib$rnaseq_depth), 2)),
               fill = "#fdcb25") +
    ylab("Mean coverage per cell") + xlab("Rank sorted position") +
    ggtitle("Mean coverage of top 5000 bases") +
    theme_classic() +
    theme(aspect.ratio = 2)
)
dev.off()



