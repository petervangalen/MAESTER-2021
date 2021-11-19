# Peter van Galen, 211108
# Visualize indel from in glioblastoma cells

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
setwd("~/DropboxMGB/Projects/Maester/AnalysisPeter/8_Glioblastoma")

# Functions & colors (available at https://github.com/vangalenlab/MAESTER-2021)
source("../210215_FunctionsGeneral.R")
heat.col <- read_excel("../MAESTER_colors.xlsx", sheet = 2, col_names = "heatcol")$heatcol
celltype.col <- c(AC = "#ffebcd", MES1 = "#ffd700", MES2 = "#bc8f8f", NPC1 = "#daa520", NPC2 = "#ff4500",
                  OPC = "#7fff00", Oligo = "#4682b4", Tumor = "#00ff7f", Malignant = "#b0e0e6", Cycling = "#5f9ea0",
                  Stromal = "#ff00ff", Myeloid = "#a9a9a9", Neutrophils = "#66cdaa", Mast = "#000000", pDCs = "#00bfff",
                  Tcells = "#696969", Bcells = "#008080", Plasma = "#9932cc")
subset.col <- c(A = "#ff69b4", C = "#4b0082", PBMC = "#66cdaa", NonMalignant = "#d3d3d3")

# Load and check data generated in 8.1_MGH252_voi.R
metadata.tib <- read_tsv(file = "MGH252_metadata.txt")
metadata.tib %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = CellType)) +
    geom_point_rast(size = 0.2) +
    scale_color_manual(values = celltype.col) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    theme_bw() +
    theme(aspect.ratio = 1, panel.grid = element_blank())

# Load Caleb's indel calls & equalize cell barcodes to metadata
A_mut.tib <- read_excel("Caleb email 211018/mut13984-glio.xlsx", sheet = 1)
C_mut.tib <- read_excel("Caleb email 211018/mut13984-glio.xlsx", sheet = 2)

A_mut.tib <- A_mut.tib %>% mutate(cell_barcode = gsub("1", "A", cell_barcode)) %>% select(-score)
C_mut.tib <- C_mut.tib %>% mutate(cell_barcode = gsub("1", "C", cell_barcode))
mut.tib <- rbind(A_mut.tib, C_mut.tib)
mut.tib <- mut.tib %>% mutate(VAF = alt_reads / (alt_reads+ref_reads)) %>% rename(cell_barcode = "cell")
mut.tib$cell %in% metadata.tib$cell

# Plot
pdf("8.2_Indel.pdf")

metadata.tib %>% left_join(mut.tib) %>% mutate(VAF = replace_na(VAF, 0)) %>% arrange(VAF) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = VAF)) +
    geom_point_rast(size = 0.5) +
    scale_color_gradientn(colors = heat.col[4:10]) +
    theme_classic() +
    theme(aspect.ratio = 1, axis.line = element_blank(), plot.title = element_text(hjust=0.5),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
    
dev.off()

