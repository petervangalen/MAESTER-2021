# Peter van Galen, 210129
# TET2 mutation frequencies in subclones. Unlike last time, I need the number of transcripts, not the number of cells.

# Prerequisites
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(tidyverse)
library(data.table)
library(Seurat)
library(Matrix)
library(gdata)

rm(list=ls())
setwd("~/DropboxPartners/Projects/Maester/AnalysisPeter/210129_TET2_mutations")

# Functions & colors
source("../201007_FunctionsGeneral.R")
popcol.df <- read.xls("~/DropboxPartners/Pipelines/AuxiliaryFiles/PopCol.xlsx", sheet = 3, row.names = 1)
mycol.ch <- popcol.df$hex
names(mycol.ch) <- rownames(popcol.df)

# Load Seurat object, extract metadata
seu <- readRDS("../210123_BPDCN712_Diagnosis/BPDCN712_Seurat.rds")

# TET2 mutated / wild-type transcript calls
xvseq.ls <- list(`TET2.S792*` = read_tsv("~/DropboxPartners/Projects/Single-cell_BPDCN/AnalysisDaniel/PCRanalysis/200108_BPDCN712/TET2.2340/TET2.2340.FilteredCells.txt"), `TET2.Q1034*` =  read_tsv("~/DropboxPartners/Projects/Single-cell_BPDCN/AnalysisDaniel/PCRanalysis/200108_BPDCN712/TET2.3078/TET2.3078.FilteredCells.txt"), `TET2.R1216*` = read_tsv("~/DropboxPartners/Projects/Single-cell_BPDCN/AnalysisDaniel/PCRanalysis/200108_BPDCN712/TET2.3626/TET2.3626.FilteredCells.txt"), `TET2.H1380Y` = read_tsv("~/DropboxPartners/Projects/Single-cell_BPDCN/AnalysisDaniel/PCRanalysis/200108_BPDCN712/TET2.4104/TET2.4104.FilteredCells.txt"))
xvseq.tib <- bind_rows(xvseq.ls, .id = "mutation") %>% mutate(cell = str_c(BC, "-1"), .before = 1) %>% select(-BC) %>%
    mutate(mutation = factor(mutation, levels = c("TET2.S792*", "TET2.Q1034*", "TET2.R1216*", "TET2.H1380Y"))) %>%
    filter(cell %in% colnames(seu))

# Add columns of which clone the cells belong to
positive_cells.tib <- read_tsv("../210123_BPDCN712_Diagnosis/210204_positive_cells.txt")
xvseq.tib <- xvseq.tib %>% left_join(select(positive_cells.tib, cell, clone), by = "cell")
write_tsv(xvseq.tib, file = "xvseq.txt")

# Summarize by clone and mutation
clones.tib <- xvseq.tib %>% group_by(clone, mutation) %>% summarize(wildtype = sum(wtUMIs), mutant = sum(mutUMIs)) %>%
    mutate(total = wildtype+mutant, fraction = mutant/(wildtype+mutant)) %>%
    mutate(clone = factor(clone, levels = unique(positive_cells.tib$clone))) %>% arrange(clone, mutation)

# Plot mutation heatmap
pdf("mutation_heatmap.pdf", height = 4, width = 12)
clones.tib %>%
    ggplot(aes(x = clone, y = mutation, fill = fraction, label = total)) +
    geom_tile() +
    geom_text() + 
    scale_y_discrete(limits = rev(c("TET2.S792*", "TET2.Q1034*", "TET2.R1216*", "TET2.H1380Y"))) +
    theme_classic() +
    theme(axis.line=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
          aspect.ratio = length(unique(clones.tib$mutation))/length(unique(clones.tib$clone)),
          plot.title = element_text(hjust = 0.5),
          axis.ticks=element_blank(), panel.background = element_rect(colour = "black", size=1, fill="#c0c0c0")) +
    scale_fill_gradient(low = "white", high = "#4B0092")
dev.off()

# Save for Tyler
clones.tib %>% pivot_wider(id_cols = mutation, names_from = clone, values_from = wildtype) %>% replace(is.na(.), 0) %>%
    write_tsv(file = "wildtype_transcripts.txt")
clones.tib %>% pivot_wider(id_cols = mutation, names_from = clone, values_from = mutant) %>% replace(is.na(.), 0) %>%
    write_tsv(file = "mutated_transcripts.txt")

# Combine all not-myeloid-biased clones
combined_clones.tib <- clones.tib %>% na.omit %>%
    mutate(clone_summary = case_when(clone == "2593_G>A" ~ "2593_G>A",
                                     clone == "6243_G>A" ~ "6243_G>A",
                                     TRUE ~ "Other_variants")) %>%
    group_by(clone_summary, mutation) %>%
    summarize(wildtype = sum(wildtype), mutant = sum(mutant), total = sum(total),
              fraction = sum(mutant)/(sum(wildtype)+sum(mutant)))

pdf("mutation_heatmap2.pdf", width = 4, height = 4)
combined_clones.tib %>%
    ggplot(aes(x = clone_summary, y = mutation, fill = fraction, label = total)) +
    geom_tile() +
    geom_text() + 
    scale_y_discrete(limits = rev(c("TET2.S792*", "TET2.Q1034*", "TET2.R1216*", "TET2.H1380Y"))) +
    theme_classic() +
    theme(axis.line=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
          aspect.ratio = 4/3,
          plot.title = element_text(hjust = 0.5),
          axis.ticks=element_blank(), panel.background =  element_rect(colour = "black", size=1, fill="#c0c0c0")) +
    scale_fill_gradient(low = "white", high = "#4B0092")
dev.off()


