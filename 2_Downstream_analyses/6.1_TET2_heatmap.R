# Peter van Galen, 210227
# TET2 mutation frequencies in subclones.

# Prerequisites
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(tidyverse)
library(data.table)
library(Seurat)
library(Matrix)
library(readxl)

rm(list=ls())
setwd("~/DropboxPartners/Projects/Maester/AnalysisPeter/6_TET2_mutations")

# Functions and colors (available at https://github.com/vangalenlab/MAESTER-2021)
source("../210215_FunctionsGeneral.R")
popcol.df <- read_excel("../MAESTER_colors.xlsx")
mycol.ch <- popcol.df$hex
names(mycol.ch) <- rownames(popcol.df)

# Load Seurat object (available at https://vangalenlab.bwh.harvard.edu/maester-2021/)
seu <- readRDS("../4_CH_sample/BPDCN712_Seurat.rds")

# TET2 mutated / wild-type transcript calls, available on https://vangalenlab.bwh.harvard.edu/maester-2021/
xvseq.ls <- list(`TET2.S792X` = read_tsv("~/DropboxPartners/Projects/Single-cell_BPDCN/AnalysisDaniel/PCRanalysis/200108_BPDCN712/TET2.2340/TET2.2340.FilteredCells.txt"),
                 `TET2.Q1034X` =  read_tsv("~/DropboxPartners/Projects/Single-cell_BPDCN/AnalysisDaniel/PCRanalysis/200108_BPDCN712/TET2.3078/TET2.3078.FilteredCells.txt"),
                 `TET2.R1216X` = read_tsv("~/DropboxPartners/Projects/Single-cell_BPDCN/AnalysisDaniel/PCRanalysis/200108_BPDCN712/TET2.3626/TET2.3626.FilteredCells.txt"),
                 `TET2.H1380Y` = read_tsv("~/DropboxPartners/Projects/Single-cell_BPDCN/AnalysisDaniel/PCRanalysis/200108_BPDCN712/TET2.4104/TET2.4104.FilteredCells.txt"))
xvseq.tib <- bind_rows(xvseq.ls, .id = "mutation") %>% mutate(cell = str_c(BC, "-1"), .before = 1) %>% select(-BC) %>%
    mutate(mutation = factor(mutation, levels = c("TET2.S792X", "TET2.Q1034X", "TET2.R1216X", "TET2.H1380Y"))) %>%
    filter(cell %in% colnames(seu))

# Add columns of which clone the cells belong to (available at https://vangalenlab.bwh.harvard.edu/maester-2021/)
positive_cells.tib <- read_tsv("../4_CH_sample/210204_positive_cells.txt")
xvseq.tib <- xvseq.tib %>% left_join(select(positive_cells.tib, cell, clone), by = "cell")
write_tsv(xvseq.tib, file = "xvseq.txt")

# Summarize by clone and mutation
clones.tib <- xvseq.tib %>% group_by(clone, mutation) %>% summarize(cell_number = n_distinct(cell), wildtype = sum(wtUMIs), mutant = sum(mutUMIs)) %>%
    mutate(total_transcripts = wildtype+mutant, fraction = mutant/(wildtype+mutant)) %>%
    mutate(clone = factor(clone, levels = unique(positive_cells.tib$clone))) %>% arrange(clone, mutation)

# Plot mutation heatmap
pdf("mutation_heatmap.pdf", height = 4, width = 12)
clones.tib %>%
    ggplot(aes(x = clone, y = mutation, fill = fraction, label = cell_number)) + # change label to total_transcripts for the number of transcripts
    geom_tile() +
    geom_text() + 
    scale_y_discrete(limits = rev(c("TET2.S792X", "TET2.Q1034X", "TET2.R1216X", "TET2.H1380Y"))) +
    theme_classic() +
    theme(axis.line=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
          aspect.ratio = length(unique(clones.tib$mutation))/length(unique(clones.tib$clone)),
          plot.title = element_text(hjust = 0.5),
          axis.ticks=element_blank(), panel.background = element_rect(colour = "black", size=1, fill="#c0c0c0")) +
    scale_fill_gradient(low = "white", high = "#4B0092")
dev.off()

# Combine all not-myeloid-biased clones
combined_clones.tib <- xvseq.tib %>% na.omit %>%
    mutate(clone_summary = case_when(clone == "2593_G>A" ~ "2593_G>A",
                                     clone == "6243_G>A" ~ "6243_G>A",
                                     TRUE ~ "Other_variants")) %>%
    group_by(clone_summary, mutation) %>%
    summarize(cell_number = n_distinct(cell), wildtype = sum(wtUMIs), mutant = sum(mutUMIs)) %>%
    mutate(total_transcripts = wildtype+mutant, fraction = mutant/(wildtype+mutant))

pdf("mutation_heatmap2.pdf", width = 4, height = 4)
combined_clones.tib %>%
    ggplot(aes(x = clone_summary, y = mutation, fill = fraction, label = cell_number)) + # change label to total_transcripts for the number of transcripts
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

#17+10+10+21+2+3+3+15+11+10+21=123 transcripts
#16+10+10+17+2+2+2+12+7+10+16=104 cells (of which 70 are unique)



# Double check if cell numbers correspond between the BPDCN project and this figure
metadata.tib <- as_tibble(seu@meta.data, rownames = "cell") %>% select(cell, contains("TET2"))
metadata.tib %>% filter(cell %in% filter(positive_cells.tib, clone == "2593_G>A")$cell) %>%
    select(-cell) %>% apply(2, function(x) sum(grepl("mutant", x), grepl("wildtype", x)))
metadata.tib %>% filter(cell %in% filter(positive_cells.tib, clone == "6243_G>A")$cell) %>%
    select(-cell) %>% apply(2, function(x) sum(grepl("mutant", x), grepl("wildtype", x)))
metadata.tib %>% filter(cell %in% filter(positive_cells.tib, ! clone %in% c("2593_G>A", "6243_G>A"))$cell) %>%
    select(-cell) %>% apply(2, function(x) sum(grepl("mutant", x), grepl("wildtype", x)))
# Perfect.
