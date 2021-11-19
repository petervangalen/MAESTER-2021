# Peter van Galen, 211107
# TET2 mutation frequencies in mtDNA subclones

# Prerequisites
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(tidyverse)
library(Seurat)
library(ComplexHeatmap)
library(circlize)

rm(list=ls())
setwd("~/DropboxMGB/Projects/Maester/AnalysisPeter/6_TET2_mutations")

# Functions and colors (available at https://github.com/vangalenlab/MAESTER-2021)
source("../210215_FunctionsGeneral.R")

# Load Seurat object (available at https://vangalenlab.bwh.harvard.edu/maester-2021/)
seu <- readRDS("../4_CH_sample/BPDCN712_Seurat_Final.rds")


# Intersect GoT with MAESTER ----------------------------------------------------------------------

# GoT mutated / wild-type transcript calls, available at https://github.com/vangalenlab/MAESTER-2021
got.ls <- list(`ASXL1.G642fs.1` = read_tsv("MutationData/ASXL1.1886.FilteredCells.txt"),
               `ASXL1.G642fs.2` = read_tsv("MutationData/ASXL1.1898.FilteredCells.txt"),
               `TET2.S792X` = read_tsv("MutationData/TET2.2340.FilteredCells.txt"),
               `TET2.Q1034X` =  read_tsv("MutationData/TET2.3078.FilteredCells.txt"),
               `TET2.R1216X` = read_tsv("MutationData/TET2.3626.FilteredCells.txt"),
               `TET2.H1380Y` = read_tsv("MutationData/TET2.4104.FilteredCells.txt"))
got.tib <- bind_rows(got.ls, .id = "mutation") %>% mutate(cell = str_c(BC, "-1"), .before = 1) %>% select(-BC) %>%
    mutate(mutation = factor(mutation, levels = c("ASXL1.G642fs.1", "ASXL1.G642fs.2", "TET2.S792X", "TET2.Q1034X",
    "TET2.R1216X", "TET2.H1380Y"))) %>% filter(cell %in% colnames(seu))
# ASXL1 was enriched with two different primers. Enumerate for Reviewer response, then collapse.
got.tib %>% filter(grepl("ASXL1", mutation), wtUMIs > 0) %>% # wtUMIs
    mutate(cell_index = as.numeric(factor(cell)), .before = 3) %>% arrange(cell_index)
got.tib <- got.tib %>% mutate(mutation = gsub(".1$|.2$", "", mutation)) %>% group_by(cell, mutation) %>%
    summarize(cell = cell[1], mutation = mutation[1], wtUMIs = max(wtUMIs), mutUMIs = max(mutUMIs)) %>% ungroup %>% unique

# Add columns of which clone the cells belong to
positive_cells.tib <- read_tsv("../4_CH_sample/4.4_positive_cells.txt")
positive_cells.tib <- positive_cells.tib %>%
    mutate(clone_summary = case_when(clone == "2593_G>A" ~ "2593_G>A",
                                     clone == "6243_G>A" ~ "6243_G>A",
                                     TRUE ~ "Other_variants"))
got.tib <- got.tib %>% left_join(select(positive_cells.tib, cell, clone_summary), by = "cell") %>% unique
write_tsv(got.tib, file = "got.txt")


# Supplemental Figure 16B -------------------------------------------------------------------------

# Denominator for top legend
positive_cells.tib %>% group_by(clone_summary) %>% summarize(n_cells = length(unique(cell)))

# Call TET2 wild-type / mutation cells and count
summary.tib <- got.tib %>% filter(mutation != "ASXL1.G642fs") %>%
    mutate(TET2.call = case_when(wtUMIs > 0 & mutUMIs == 0 ~ "wt", mutUMIs > 0 ~ "mut")) %>%
    group_by(clone_summary, mutation) %>% summarize(wt_cells = sum(TET2.call == "wt"),
                                                    mut_cells = sum(TET2.call == "mut"),
                                                    total_cells = n()) %>%
    na.omit()

# Numerator for top legend
summary.tib %>% group_by(clone_summary) %>% summarize(total_cells = sum(total_cells))
# Cell number for left legend
summary.tib %>% group_by(mutation) %>% summarize(total_cells = sum(total_cells))

pdf("mutation_heatmap.pdf", width = 6, height = 6)
summary.tib %>%
    ggplot(aes(x = clone_summary, y = mutation, fill = mut_cells / total_cells,
               label = str_c(wt_cells, " wt\n", mut_cells, " mut"))) +
    geom_tile() +
    geom_text() + 
    scale_y_discrete(limits = rev(c("TET2.S792X", "TET2.Q1034X", "TET2.R1216X", "TET2.H1380Y"))) +
    theme_classic() +
    theme(axis.line=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
          aspect.ratio = 4/3,
          plot.title = element_text(hjust = 0.5),
          axis.ticks=element_blank(), panel.background =  element_rect(colour = "black", size=1, fill="#c0c0c0")) +
    scale_fill_gradient(low = "white", high = "#4B0092")
dev.off()


# Supplemental Figure 16C -------------------------------------------------------------------------

# Number of positive cells for each mutation
mut.ls <- list(X2593GA = positive_cells.tib %>% filter(clone == "2593_G>A") %>% .$cell %>% unique,
               ASXL1.G642fs = got.tib %>% filter(mutation == "ASXL1.G642fs" & mutUMIs > 0) %>% .$cell %>% unique,
               TET2.S792X = got.tib %>% filter(mutation == "TET2.S792X" & mutUMIs > 0) %>% .$cell %>% unique,
               TET2.Q1034X = got.tib %>% filter(mutation == "TET2.Q1034X" & mutUMIs > 0) %>% .$cell %>% unique,
               TET2.R1216X = got.tib %>% filter(mutation == "TET2.R1216X" & mutUMIs > 0) %>% .$cell %>% unique,
               TET2.H1380Y = got.tib %>% filter(mutation == "TET2.H1380Y" & mutUMIs > 0) %>% .$cell %>% unique)
lengths(mut.ls)

# Intersect all
mut.ls.mat <- as.matrix(do.call(cbind, lapply(mut.ls, function(x) lapply(mut.ls, function(y) length(intersect(x, y))))))
mut.mat <- apply(mut.ls.mat, 1, function(x) as.numeric(x))
rownames(mut.mat) <- rownames(mut.ls.mat)

col_fun <- colorRamp2(c(0, 1, max(mut.mat)), c("#DDDDDD", "white", "#4B0092"))

pdf("cooccurrence.pdf", width = 5, height = 4)
Heatmap(mut.mat,
        cluster_rows = F,
        cluster_columns = F,
        col = col_fun,
        rect_gp = gpar(col = "black", lwd = 1),
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(mut.mat[i, j], x, y, gp = gpar(fontsize = 12))
        })
dev.off()






