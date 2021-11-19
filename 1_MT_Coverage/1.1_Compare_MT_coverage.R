# Peter van Galen, 210412
# Visualize chrM coverage from maegatk data objects. This script plots two maegatk objects together. For just one object, see 1.2_MT_coverage.R instead.

options(scipen = 999)

library(tidyverse)
library(Matrix)
library(ggforce)
library(SummarizedExperiment)

rm(list=ls())
setwd("~/DropboxMGB/Projects/Maester/AnalysisPeter/1_MT_Coverage")

# Functions (available at https://github.com/vangalenlab/MAESTER-2021)
source("../210215_FunctionsGeneral.R")


#~~~~~~~~~~~#
# Load data #
#~~~~~~~~~~~#

# Locally saved maegatk data (available at https://vangalenlab.bwh.harvard.edu/maester-2021/) - choose one
    experiment.name <- "SW_CellLineMix"; method1 <- "RNAseq alone"; method2 <- "RNAseq+MAESTER"
    method1.maegatk <- readRDS(file = paste0(experiment.name, "_RNAseq_mr3_maegatk.rds"))
    method2.maegatk <- readRDS(file = paste0(experiment.name, "_All_mr3_maegatk.rds"))

    experiment.name <- "TenX_CellLineMix"; method1 <- "RNAseq alone"; method2 <- "RNAseq+MAESTER"
    method1.maegatk <- readRDS(file = paste0(experiment.name, "_RNAseq_mr3_maegatk.rds"))
    method2.maegatk <- readRDS(file = paste0(experiment.name, "_All_mr3_maegatk.rds"))

    experiment.name <- "TenX_BPDCN712"; method1 <- "RNAseq alone"; method2 <- "RNAseq+MAESTER"
    method1.maegatk <- readRDS(file = paste0(experiment.name, "_RNAseq_mr3_maegatk.rds"))
    method2.maegatk <- readRDS(file = paste0(experiment.name, "_All_mr3_maegatk.rds"))

    experiment.name <- "TenX_Multiome_LLL"; method1 <- "MAESTER"; method2 <- "ATAC-LLL"
    method1.maegatk <- readRDS(file = "TenX_Multiome_LLL_mr3_maegatk.rds")
    method2.maegatk <- readRDS(file = "PBMC_LLL_atac_mgatk.rds")

#    experiment.name <- "TenX_Multiome_DIG"; method1 <- "ATAC-DIG"; method2 <- "MAESTER"
#    method1.maegatk <- readRDS(file = "PBMC_DIG_atac_mgatk.rds")
#    method2.maegatk <- readRDS(file = "TenX_Multiome_DIG_mr3_maegatk.rds")


#~~~~~~#
# Plot #
#~~~~~~#

pdf(file = paste0("1.1_", experiment.name, ".pdf"))

# Plot coverage per position ----------------------------------------------------------------------

# Set y axis parameters
if (experiment.name %in% c("SW_CellLineMix", "TenX_Multiome_DIG", "TenX_Multiome_LLL")) {
    ymax <- 200
} else if (experiment.name %in% c("TenX_CellLineMix", "TenX_BPDCN712")) {
    ymax <- 1000
}

# Gene locations
GenePos.tib <- tibble(Names = c("MT.ATP6", "MT.ATP8", "MT.CO1", "MT.CO2", "MT.CO3", "MT.CYB", "MT.ND1", "MT.ND2", "MT.ND3",
                                "MT.ND4", "MT.ND4L", "MT.ND5", "MT.ND6", "MT.RNR1", "MT.RNR2"),
                      start = c(8527, 8366, 5904, 7586, 9207, 14747, 3307, 4470, 10059, 10760, 10470, 12337, 14149, 648, 1671), 
                      end = c(9207, 8572, 7445, 8269, 9990, 15887, 4262, 5511, 10404, 12137, 10766, 14148, 14673, 1601, 3229))
GenePos.tib <- GenePos.tib %>% arrange(start) %>%
    mutate(mid = round((end-start)/2+start,0), ycoord = rep(c(ymax*1.2,ymax*1.1), length.out = 15))

# Plot
base.tib <- tibble(base = 1:16569,
                   method1_depth = rowMeans(assays(method1.maegatk)[["coverage"]]),
                   method2_depth = rowMeans(assays(method2.maegatk)[["coverage"]]))


print(
    base.tib %>% ggplot() +
        geom_bar(aes(x = base, y = ifelse(method2_depth > 1, yes = method2_depth, no = NA)), stat = "identity", fill = "#64b53b", width = 1) + 
        geom_bar(aes(x = base, y = ifelse(method1_depth > 1, yes = method1_depth, no = NA)), stat = "identity", fill = "#fdcb25", width = 1) +
        coord_cartesian(ylim = c(1, ymax), xlim = c(700, 15900)) +
        scale_y_continuous(trans = "log10") +
        geom_segment(data = GenePos.tib, aes(x = start, y = ycoord, xend = end, yend = ycoord)) +
        geom_text(data = GenePos.tib, aes(x = mid, y = ycoord-ymax*0.2, label = cutf(Names, d = "\\.", f = 2)), size = 3) +
        ylab("Mean coverage per cell") + xlab("Position along chrM") +
        theme_classic() +
        theme(aspect.ratio = 0.5)
)

# Plot mean coverage for top 500 cells ------------------------------------------------------------

# Use same cells for coverage from RNA-seq and Maester
common.cells <- intersect(colnames(method1.maegatk), colnames(method2.maegatk))
method1.maegatk <- method1.maegatk[,common.cells]
method2.maegatk <- method2.maegatk[,common.cells]

# Mean coverage per cell
cells.tib <- tibble(cell = colnames(method1.maegatk),
                    method1_depth = method1.maegatk$depth,
                    method2_depth = method2.maegatk$depth)
topcells.tib <- cells.tib %>% slice_max(order_by = method1_depth + method2_depth, n = 500)

print(
pivot_longer(topcells.tib, cols = c("method1_depth", "method2_depth")) %>%
    mutate(name = gsub("method1_depth", method1, gsub("method2_depth", method2, name))) %>%
    mutate(name = factor(name, levels = c(method1, method2))) %>%
    ggplot(aes(x = name, y = value, color = name)) +
    geom_sina(size = 0.3) +
    scale_color_manual(values = c("#fdcb25", "#64b53b")) +
    geom_violin(color = "black", fill = NA) +
    coord_cartesian(ylim = c(0.1, 800)) +
    scale_y_continuous(trans = "log10") +
    ylab("Mean coverage per cell") + xlab("") +
    annotate("text", x = c(1,2), y = c(max(topcells.tib$method1_depth),max(topcells.tib$method2_depth))*1.5,
             label = round(c(mean(topcells.tib$method1_depth), mean(topcells.tib$method2_depth)), 2)) +
    theme_classic() +
    theme(aspect.ratio = 2, plot.title = element_text(hjust = 0.5), legend.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    ggtitle("Mean coverage of top 500 cells")
)

# Plot mean depth for top 5000 bases --------------------------------------------------------------
top.tib <- base.tib %>% arrange(desc(method2_depth)) %>% mutate(key = row_number(), .before = 1) %>% filter(key %in% 1:5000)
print(
ggplot(top.tib) +
    geom_bar(aes(x = key, y = method2_depth), stat = "identity", fill = "#64b53b", width = 1) +
    geom_bar(aes(x = key, y = ifelse(method1_depth > 1, yes = method1_depth, no = NA)), stat = "identity", fill = "#fdcb25", width = 1) +
    coord_cartesian(ylim = c(1, ymax)) +
    scale_y_continuous(trans = "log10") +
    geom_label(data = data.frame(), aes(x = 2500, y = mean(top.tib$method2_depth), label = round(mean(top.tib$method2_depth), 2)),
               fill = "#64b53b") +
    geom_label(data = data.frame(), aes(x = 2500, y = 2, label = round(mean(top.tib$method1_depth), 2)),
               fill = "#fdcb25") +
    ylab("Mean coverage per cell") + xlab("Rank sorted position") +
    ggtitle("Mean coverage of top 5000 bases") +
    theme_classic() +
    theme(aspect.ratio = 2)
)
dev.off()



