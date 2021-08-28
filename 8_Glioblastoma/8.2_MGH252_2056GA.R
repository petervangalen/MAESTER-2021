# Peter van Galen, 210823
# Identify informative variants in glioblastoma sample

#~~~~~~~~~~~~~~~#
# Prerequisites #
#~~~~~~~~~~~~~~~#

options(stringsAsFactors = FALSE)
options(scipen = 999)

library(tidyverse)
library(data.table)
library(SummarizedExperiment)
library(readxl)
library(ggrastr)
library(Matrix)
library(GGally) # for ggpairs
library(ggrepel)
library(gridExtra)

rm(list=ls())
setwd("~/DropboxPartners/Projects/Maester/AnalysisPeter/8_Glioblastoma")

## Functions & colors (available at https://github.com/vangalenlab/MAESTER-2021)
source("../210215_FunctionsGeneral.R")
heat.col <- read_excel("../MAESTER_colors.xlsx", sheet = 2, col_names = "heatcol")$heatcol
celltype.col <- c(AC = "#ffebcd", MES1 = "#ffd700", MES2 = "#bc8f8f", NPC1 = "#daa520", NPC2 = "#ff4500",
                  OPC = "#7fff00", Oligo = "#4682b4", Tumor = "#00ff7f", Malignant = "#b0e0e6", Cycling = "#5f9ea0",
                  Stromal = "#ff00ff", Myeloid = "#a9a9a9", Neutrophils = "#66cdaa", Mast = "#000000", pDCs = "#00bfff",
                  Tcells = "#696969", Bcells = "#008080", Plasma = "#9932cc")
subset.col <- c(A = "#ff69b4", C = "#4b0082", PBMC = "#66cdaa", NonMalignant = "#d3d3d3")


#~~~~~~~~~~~~~~~~~~~~~#
# Load + wrangle data #
#~~~~~~~~~~~~~~~~~~~~~#

# Wrangle Chadi's metadata  (available at https://github.com/vangalenlab/MAESTER-2021)
metadata.full <- read_tsv("Chadi Slack 210819/MGH252_NoM_A_C_PBMC_MetaData.txt")
umap.full <- read_tsv("Chadi Slack 210819/MGH252_NoM_A_C_PBMC_UMAP.txt")
metadata.full <- metadata.full %>% left_join(umap.full) %>%
    mutate(Site = gsub("[0-9]", "", cutf(X1, d = "\\.|_", f = 4)), .before = 2) %>%
    mutate(cell = str_c(cutf(X1, d = "\\.", f = 2), "-", Site), .before = 2) %>%
    select(cell, Site, CellType, UMAP_1, UMAP_2)

# Wrangle maegatk data (generated in script 8.1_MGH252_voi.R)
maegatk <- readRDS(file = "MGH252_maegatk.rds")
metadata.tib <- filter(metadata.full, cell %in% colnames(maegatk))
maegatk <- maegatk[,metadata.tib$cell]
af.dm <- data.matrix(computeAFMutMatrix(maegatk))*100
all(colnames(af.dm) == metadata.tib$cell)

# Add variant info. 2056_G>A is for the paper supplement
var <- "2056_G>A"
metadata.tib$vaf <- af.dm[var,]
metadata.tib$cov <- assays(maegatk)[["coverage"]][as.numeric(cutf(var, d = "_")),]


#~~~~~~~~~~~#
# Visualize #
#~~~~~~~~~~~#

g1 <- metadata.tib %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = CellType)) +
    geom_point_rast(size = 0.5) +
    scale_color_manual(values = celltype.col) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    theme_bw() +
    theme(aspect.ratio = 1, panel.grid = element_blank())

g2 <- metadata.tib %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = Site)) +
    geom_point_rast(size = 0.5) +
    scale_color_manual(values = subset.col) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    theme_bw() +
    theme(aspect.ratio = 1, panel.grid = element_blank())

g3 <- metadata.tib %>%
    filter(cov >= 5) %>%
    arrange(vaf) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = vaf)) +
    geom_point_rast(size = 0.5) +
    scale_color_gradientn(colors = heat.col[4:10], n.breaks = 3) +
    theme_classic() +
    theme(aspect.ratio = 1, axis.line = element_blank(), plot.title = element_text(hjust=0.5),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    ggtitle(var)

pdf("8.2_1_Plots_2056.pdf", width = 15, height = 6)
grid.arrange(g1, g2, g3, nrow = 1)
dev.off()

