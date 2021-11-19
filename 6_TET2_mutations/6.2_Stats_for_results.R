# Peter van Galen, 211107
# Various numbers for results text and response to reviewers. Much of this is superfluous.

# Prerequisites
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(tidyverse)
library(Seurat)
library(SummarizedExperiment)
library(Matrix)
library(readxl)

rm(list=ls())
setwd("~/DropboxMGB/Projects/Maester/AnalysisPeter/6_TET2_mutations")

# Functions and colors (available at https://github.com/vangalenlab/MAESTER-2021)
source("../210215_FunctionsGeneral.R")
popcol.df <- read_excel("../MAESTER_colors.xlsx")
mycol.ch <- popcol.df$hex
names(mycol.ch) <- rownames(popcol.df)

# Load Seurat object (available at https://vangalenlab.bwh.harvard.edu/maester-2021/)
seu <- readRDS("../4_CH_sample/BPDCN712_Seurat_Final.rds")

# Load Maegtk, calculate allele frequencies (available at https://vangalenlab.bwh.harvard.edu/maester-2021/)
maegatk.rse <- readRDS("../4_CH_sample/BPDCN712_Maegatk_Final.rds")
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))*100

# Clone information (available at https://github.com/vangalenlab/MAESTER-2021)
positive_cells.tib <- read_tsv("../4_CH_sample/4.4_positive_cells.txt")

# TET2 and ASXL1 mutated / wild-type transcript calls (generated in the script 6.1_TET2_heatmap.R)
got.tib <- read_tsv("got.txt")


# Get stats for the results section ---------------------------------------------------------------

ncol(seu)
# 9346 cells in total

### GoT - All cells ###
sum(got.tib$wtUMIs)
# 509 wild-type transcripts in total

sum(got.tib$mutUMIs)
# 55 mutated transcripts in total

got.tib %>% filter(wtUMIs > 0) %>% .$cell %>% unique %>% length
# 347 cells with one or more wild-type transcripts

got.tib %>% filter(mutUMIs > 0) %>% .$cell %>% unique %>% length
# 45 cells with one or more mutated transcripts

got.tib$cell %>% unique %>% length
# 374 cells with any genotyping information. 374/9346*100 = 4.00

got.tib %>% filter(grepl("ASXL1", mutation)) %>% .$cell %>% unique %>% length
got.tib %>% filter(grepl("TET2", mutation)) %>% .$cell %>% unique %>% length
# 39 with ASXL1, 338 for TET2 (39/9346*100 = 0.4%, 338/9346*100 = 3.6%)

### GoT - Subclones ###
positive_cells.tib %>% .$cell %>% unique %>% length
# 1,397 cells that are part of any clone

# Filter got data by these cells
got_clones.tib <- got.tib %>% filter(cell %in% positive_cells.tib$cell)

sum(got_clones.tib$wtUMIs)
# 118 wild-type transcripts in total

sum(got_clones.tib$mutUMIs)
# 17 mutated transcripts

got_clones.tib %>% filter(wtUMIs > 0) %>% .$cell %>% unique %>% length
# 73 cells with one or more wild-type transcripts

got_clones.tib %>% filter(mutUMIs > 0) %>% .$cell %>% unique %>% length
# 15 cells with one or more mutated transcripts

got_clones.tib$cell %>% unique %>% length
# 81 cells with any genotyping information. 
got_clones.tib %>% filter(grepl("TET2", mutation)) %>% .$cell %>% unique %>% length
# 70 cells with TET2 genotyping information. 70/1397*100 = 5.01%. The following yields the same result:
unique(positive_cells.tib$cell) %in% filter(got.tib, grepl("TET2", mutation))$cell %>% mean * 100 # 5.01% genotyping within clones
colnames(seu) %in% filter(got.tib, grepl("TET2", mutation))$cell %>% mean * 100 # 3.6% overall
# This is less than the sum of numbers in Supplemental Figure 13B because some cells have multiple transcripts and are thus represented in multiple rows.
# ASXL1 genotpying efficiency:
colnames(seu) %in% filter(got.tib, grepl("ASXL1", mutation))$cell %>% mean * 100 # 0.4% overall
unique(positive_cells.tib$cell) %in% filter(got.tib, grepl("ASXL1", mutation))$cell %>% mean * 100 # 0.9% within clones

transcripts_per_cell.num <- apply(got_clones.tib, 1, function(x) sum(as.numeric(x[3]), as.numeric(x[4])))
names(transcripts_per_cell.num) <- got_clones.tib$cell
sum(transcripts_per_cell.num) # 135 transcripts
length(unique(names(transcripts_per_cell.num))) # 81 cells

### MAESTER - All cells ###
# Cell IDs, assays
maegatk.rse$sample %>% length

# Extract coverage from the matrix with position x cells
coverage_2593 <- assays(maegatk.rse)[["coverage"]][2593,]
par(pty="s")
plot(sort(coverage_2593), log = "y", xlab = "rank sorted cells", ylab = "2593 coverage")
length(coverage_2593); sum(coverage_2593 > 0); 9339/9346*100

tib.2593 <- tibble(cell = colnames(maegatk.rse),
                   coverage = coverage_2593,
                   A_counts_fw = assays(maegatk.rse)[["A_counts_fw"]][2593,],
                   A_counts_rev = assays(maegatk.rse)[["A_counts_rev"]][2593,],
                   C_counts_fw = assays(maegatk.rse)[["C_counts_fw"]][2593,],
                   C_counts_rev = assays(maegatk.rse)[["C_counts_rev"]][2593,],
                   G_counts_fw = assays(maegatk.rse)[["G_counts_fw"]][2593,],
                   G_counts_rev = assays(maegatk.rse)[["G_counts_rev"]][2593,],
                   T_counts_fw = assays(maegatk.rse)[["T_counts_fw"]][2593,],
                   T_counts_rev = assays(maegatk.rse)[["T_counts_rev"]][2593,])

tib.2593 %>% select(contains("counts")) %>% colSums
# 9 + 100107 = 100116 wt transcripts (G)
# 5203 mutant transcripts (A)

sum(coverage_2593)
tib.2593 %>% select(contains("counts")) %>% sum()
# 105325 total transcripts

tib.2593 %>% select(contains("counts")) %>% apply(., 2, function(x) sum(x > 0))
# 269 mutant cells

tib.2593 %>% filter(G_counts_fw > 0 | G_counts_rev > 0) %>% nrow
tib.2593 %>% filter(G_counts_fw == 0 & G_counts_rev == 0) %>% nrow
# 9300 wildtype cells (46 cells do not have any wild-type reads)

sum(tib.2593 %>% select(contains("counts")) %>% rowSums > 0)
sum(coverage_2593 > 0)
# 9339 cells with some coverage. 9339/9346*100 = 99.9%


### MAESTER - Subclones ###
maegatk_clones <- maegatk.rse[,unique(positive_cells.tib$cell)]

tib.clones.2593 <- tibble(cell = colnames(maegatk_clones),
                          coverage = assays(maegatk_clones)[["coverage"]][2593,],
                          A_counts_fw = assays(maegatk_clones)[["A_counts_fw"]][2593,],
                          A_counts_rev = assays(maegatk_clones)[["A_counts_rev"]][2593,],
                          C_counts_fw = assays(maegatk_clones)[["C_counts_fw"]][2593,],
                          C_counts_rev = assays(maegatk_clones)[["C_counts_rev"]][2593,],
                          G_counts_fw = assays(maegatk_clones)[["G_counts_fw"]][2593,],
                          G_counts_rev = assays(maegatk_clones)[["G_counts_rev"]][2593,],
                          T_counts_fw = assays(maegatk_clones)[["T_counts_fw"]][2593,],
                          T_counts_rev = assays(maegatk_clones)[["T_counts_rev"]][2593,])

tib.clones.2593 %>% select(contains("counts")) %>% colSums
# 1 + 16563 = 16564 wt transcripts (G)
# 5201 mutant transcripts (A)

tib.clones.2593 %>% select(contains("counts")) %>% sum()
# 21767 total transcripts

tib.clones.2593 %>% select(contains("counts")) %>% apply(., 2, function(x) sum(x > 0))
# 267 mutant cells

tib.clones.2593 %>% filter(G_counts_fw > 0 | G_counts_rev > 0) %>% nrow
tib.clones.2593 %>% filter(G_counts_fw == 0 & G_counts_rev == 0) %>% nrow
# 1357 wildtype cells (40 cells do not have any wild-type reads; These are likely cells with high heteroplasmy)

sum(tib.clones.2593 %>% select(contains("counts")) %>% rowSums > 0)
# 1396 cells with some coverage. 1396/1397*100 = 99.9%

tib.clones.2593$coverage %>% mean
# Mean coverage of 15.6


# Get more stats for Reviewer 1 -------------------------------------------------------------------

# What is the allele frequency of 2593_G>A ?
mean( af.dm["2593_G>A",] ) # The mean VAF per single cell is 1.27%
wt.num <- sum(tib.2593$G_counts_fw, tib.2593$G_counts_rev)
mut.num <- sum(tib.2593$A_counts_fw, tib.2593$A_counts_rev)
mut.num / (wt.num + mut.num) * 100 # More importantly, the overall VAF is 4.94%

# Wrangle
got2.tib <- got.tib %>% select(-wtUMIs) %>% mutate(mutUMIs = gsub(0, NA, mutUMIs)) %>%
    pivot_wider(names_from = c(mutation), values_from = mutUMIs)

# How many cells were genotyped (ASXL1, TET2)
got.tib %>% filter(grepl("TET2", mutation)) %>% .$cell %>% unique %>% length # 338; 338/9346*100=3.6%
got.tib %>% filter(grepl("ASXL1", mutation)) %>% .$cell %>% unique %>% length # 39; 39/9346*100=0.4%

# Remove all clones but 2593_G>A
got2.tib <- got2.tib %>% mutate(clone = factor(clone_summary, levels = c("2593_G>A"))) # "6243_G>A",

got2.tib$na_count <- apply(got2.tib, 1, function(x) sum(is.na(x)))
    
got2.tib %>% filter(na_count < 5) %>% .[,c("cell", "ASXL1.G642fs", "clone", "TET2.S792X", "TET2.Q1034X",
                                           "TET2.R1216X", "TET2.H1380Y")] %>%
    view

# ASXL1 stats for R1Q2 & R1Q3
asxl1_mut.cells <- got.tib %>% filter(grepl("ASXL1", mutation), mutUMIs > 0) %>% .$cell
got.tib %>% filter(cell %in% asxl1_mut.cells) %>% .$cell %>% unique %>% length
got.tib %>% filter(cell %in% asxl1_mut.cells) %>% mutate(cell_index = as.numeric(factor(cell)), .before = 1) %>%
    arrange(cell_index) %>% view
length(unique(asxl1_mut.cells))

# ASXL1 was genotyped in 0.4% of cells
asxl1_wt.cells <- got.tib %>% filter(grepl("ASXL1", mutation), wtUMIs > 0) %>% .$cell
length(unique(asxl1_wt.cells))
length(unique(c(asxl1_mut.cells, asxl1_wt.cells))) / ncol(seu) * 100



