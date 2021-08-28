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
seu <- readRDS("../4_CH_sample/BPDCN712_Seurat_Final.rds")

# TET2 mutated / wild-type transcript calls, available at https://github.com/vangalenlab/MAESTER-2021
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

# Add columns of which clone the cells belong to (available at https://vangalenlab.bwh.harvard.edu/maester-2021/)
positive_cells.tib <- read_tsv("../4_CH_sample/4.3_positive_cells.txt")
got.tib <- got.tib %>% left_join(select(positive_cells.tib, cell, clone), by = "cell") %>% unique
write_tsv(got.tib, file = "got.txt")

# Summarize by clone and mutation
clones.tib <- got.tib %>% group_by(clone, mutation) %>%
    summarize(cell_number = n_distinct(cell), wildtype = sum(wtUMIs), mutant = sum(mutUMIs)) %>%
    mutate(total_transcripts = wildtype+mutant, fraction = mutant/(wildtype+mutant)) %>%
    mutate(clone = factor(clone, levels = unique(positive_cells.tib$clone))) %>% arrange(clone, mutation)

# Plot mutation heatmap
pdf("mutation_heatmap.pdf", height = 4, width = 12)
clones_tet.tib <- clones.tib %>% filter(! grepl("ASXL1", mutation)) %>%
    mutate(mutation = factor(mutation, levels = c("TET2.S792X", "TET2.Q1034X", "TET2.R1216X", "TET2.H1380Y")))
clones_tet.tib %>%
    ggplot(aes(x = clone, y = mutation, fill = fraction, label = cell_number)) + # change label to total_transcripts for the number of transcripts
    geom_tile() +
    geom_text() + 
    scale_y_discrete(limits = rev(c("TET2.S792X", "TET2.Q1034X", "TET2.R1216X", "TET2.H1380Y"))) +
    theme_classic() +
    theme(axis.line=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
          aspect.ratio = length(unique(clones_tet.tib$mutation))/length(unique(clones_tet.tib$clone)),
          plot.title = element_text(hjust = 0.5),
          axis.ticks=element_blank(), panel.background = element_rect(colour = "black", size=1, fill="#c0c0c0")) +
    scale_fill_gradient(low = "white", high = "#4B0092")
dev.off()

# Combine all not-myeloid-biased clones
combined_clones.tib <- got.tib %>% na.omit %>%
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Overlap between genotyped cells $
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Two cells have both TET2.S792X and TET2.Q1034X, but they did not pass RNA-seq QC, so they were excluded.
intersect(got.ls[[3]][got.ls[[3]]$mutUMIs > 0,]$BC, got.ls[[4]][got.ls[[4]]$mutUMIs > 0,]$BC) %in% cutf(colnames(seu), d = "-")
# One cell that has TET2.R1216X and TET2.H1380Y will remain
intersect(got.ls[[5]][got.ls[[5]]$mutUMIs > 0,]$BC, got.ls[[6]][got.ls[[6]]$mutUMIs > 0,]$BC) %in% cutf(colnames(seu), d = "-")

# There is no overlap between ASXL1 and TET2 mutated cells
names(got.ls)
for (x in 1:2) { 
    for (y in 3:5) { print(
        intersect(got.ls[[x]][got.ls[[x]]$mutUMIs > 0,]$BC, got.ls[[y]][got.ls[[y]]$mutUMIs > 0,]$BC) %in% cutf(colnames(seu), d = "-")
) } }

