# Peter van Galen, 210204
# Assess lineage bias in Patient 10 diagnosis clones


#~~~~~~~~~~~~~~~~~~~~~#
#### Prerequisites ####
#~~~~~~~~~~~~~~~~~~~~~#

options(stringsAsFactors = FALSE)
options(scipen = 999)

library(tidyverse)
library(data.table)
library(SummarizedExperiment)
library(Seurat)
library(Matrix)
library(ggrastr)
library(ComplexHeatmap)
library(gdata)
library(readxl)
library(ggrepel)
#library(circlize)
#library(GGally) # for ggpairs
#library(ggforce) # for geom_sina
#library(ggbeeswarm)
#library(stringr)

rm(list=ls())
setwd("~/DropboxPartners/Projects/Maester/AnalysisPeter/210123_BPDCN712_Diagnosis")

# Functions & colors
source("../201007_FunctionsGeneral.R")
popcol.df <- read.xls("~/DropboxPartners/Pipelines/AuxiliaryFiles/PopCol.xlsx", sheet = 3)
cell_type_colors <- deframe(popcol.df[1:14,1:2])
clone_colors <- deframe(popcol.df[15:37,1:2])

### Import data
# Load Seurat object
seu <- readRDS("BPDCN712_Seurat.rds")

# Load Maegtk, calculate allele frequencies
maegatk.rse <- readRDS("BPDCN712_Maegatk.rds")
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))*100

all(colnames(af.dm) == colnames(seu))

# Extract Seurat metadata
cells.tib <- as_tibble(seu@meta.data, rownames = "cell")

# Variants of interest
voi.ch <- read_tsv("210124_Dx_vois.txt")$var


#~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Cell type numbers ####
#~~~~~~~~~~~~~~~~~~~~~~~~~#

# Make a list of cell ids that are positive for each of the variants of interest
positive_cells.ls <- list()

for (v in voi.ch) {
    message(v)
    # Determine cells with an appreciable VAF
    current_cells.ch <- colnames(af.dm)[af.dm[v,]>1]
    # Save cell IDs for positive cells
    positive_cells.ls[[v]] <- current_cells.ch
}

# Save a tibble with cell, variant, clone, cell type information
positive_cells.tib <- as_tibble(bind_rows(lapply(positive_cells.ls, function(x) data.frame(cell = x)), .id = "variant")[,2:1]) %>%
    mutate(variant = factor(variant, levels = voi.ch))
positive_cells.tib <- positive_cells.tib %>%
    mutate(clone = gsub("6205_G>A|9164_T>C", "6205_G>A-9164_T>C", gsub(
        "10158_T>A|6293_T>C", "10158_T>A-6293_T>C", gsub(
            "1415_G>A|9753_G>A", "1415_G>A-9753_G>A", variant)))) %>%
    mutate(clone = factor(clone, levels = unique(clone)))
positive_cells.tib <- positive_cells.tib %>% left_join(dplyr::select(as_tibble(seu@meta.data, rownames = "cell"), cell, CellType))
write_tsv(positive_cells.tib, file = "210204_positive_cells.txt")

# Save a tibble with cell type counts per clone
cell_type_numbers.tib <- cells.tib %>% group_by(CellType) %>% dplyr::count(name = "all_cells") %>% ungroup
cell_type_numbers_variants.tib <- positive_cells.tib %>% select(variant, CellType) %>% group_by(variant, CellType) %>% summarize(number = n()) %>%
    pivot_wider(id_cols = CellType, names_from = variant, values_from = number) %>% replace(is.na(.), 0)
cell_type_numbers.tib <- left_join(cell_type_numbers.tib, cell_type_numbers_variants.tib, by = "CellType")
write_tsv(cell_type_numbers.tib, file = "210204_cell_type_numbers.txt")


#~~~~~~~~~~~~~~~~#
#### Bar plot ####
#~~~~~~~~~~~~~~~~#

# Transpose, pivot
transpose_numbers.tib <- cell_type_numbers.tib %>% pivot_longer(c(all_cells, all_of(voi.ch)), names_to = "clone") %>% pivot_wider(names_from = CellType)

plot.tib <- transpose_numbers.tib %>%
    mutate(clone = factor(clone, levels = c("all_cells", setdiff(clone, "all_cells")))) %>%
    pivot_longer(-clone, names_to = "CellType", values_to = "Number of cells") %>%
    mutate(CellType = factor(CellType, levels = popcol.df$name[1:14]))

pdf("210204_1_CellTypeBarplots.pdf")

ggplot(filter(plot.tib, clone != "all_cells"), aes(x = clone, y = `Number of cells`, fill = CellType)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = cell_type_colors) +
    scale_y_continuous(limits=c(0, 300)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 0.7)

ggplot(plot.tib, aes(x = clone, y = `Number of cells`, fill = CellType)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = popcol.df[plot.tib$CellType,"hex"]) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 0.7)

dev.off()

# Add frequency and normalized frequency
cell_type_freq.tib <- cell_type_numbers.tib %>%
    mutate(across(-CellType, .fns = list(freq = ~ .x / sum(.x) * 100,           # equivalent to function(x) {x / sum(x) * 100 }
                                         norm = ~ (.x / sum(.x)) / (all_cells / sum(all_cells)))))

# Save for Tyler's P-value calcuations
cell_type_freq.tib %>% dplyr::select(CellType, contains("freq")) %>% write_tsv(file = "210204_cell_type_freq.txt")
cell_type_freq.tib %>% dplyr::select(CellType, contains("norm")) %>% write_tsv(file = "210204_cell_type_norm.txt")


#~~~~~~~~~~~~~~~~~~#
#### Radar plot ####
#~~~~~~~~~~~~~~~~~~#

# Now, use clones instead of variants
cell_type_numbers2.tib <- cells.tib %>% group_by(CellType) %>% dplyr::count(name = "all_cells") %>% ungroup
cell_type_numbers_clones.tib <- positive_cells.tib %>% filter(!duplicated(cell)) %>%
    select(clone, CellType) %>% group_by(clone, CellType) %>% summarize(number = n()) %>%
    pivot_wider(id_cols = CellType, names_from = clone, values_from = number) %>% replace(is.na(.), 0)
cell_type_numbers2.tib <- left_join(cell_type_numbers2.tib, cell_type_numbers_clones.tib, by = "CellType")

# Add frequency and normalized frequency
cell_type_freq2.tib <- cell_type_numbers2.tib %>%
    mutate(across(-CellType, .fns = list(freq = ~ .x / sum(.x) * 100,           # equivalent to function(x) {x / sum(x) * 100 }
                  norm = ~ (.x / sum(.x)) / (all_cells / sum(all_cells)))))

# Which ones to plot?
current.clones <- cell_type_freq2.tib %>% select(-CellType, -contains("freq"), -contains("norm")) %>% colSums
current.clones <- c(names(current.clones)[current.clones > 50], "6243_G>A")

# Function for coord_polar with straight lines (from https://stackoverflow.com/questions/57209060/how-to-draw-a-radar-plot-in-ggplot-using-polar-coordinates) and http://www.cmap.polytechnique.fr/~lepennec/en/post/radar/radarandparallelplots/
coord_radar <- function (theta = "x", start = 0, direction = 1) {
    theta <- match.arg(theta, c("x", "y"))
    r <- if (theta == "x") "y" else "x"
    ggproto("CoordRadar", CoordPolar, theta = theta, r = r, start = start, 
            direction = sign(direction),
            is_linear = function(coord) TRUE)
}

# Plot frequency
cell_type_freq2.tib %>% dplyr::select(CellType, contains(current.clones)) %>%
    add_row(CellType = "empty") %>%
    replace(is.na(.), 0) %>%
    mutate(CellType = factor(CellType, levels = CellType)) %>%
    dplyr::select(CellType, contains("freq")) %>% pivot_longer(cols = -CellType) %>%
    mutate(name = factor(name, levels = str_c(current.clones, "_freq"))) %>%
    ggplot(aes(x = as.numeric(CellType), y = value, fill = name, color = name)) +
    coord_radar() +
    geom_line() +
    geom_polygon(alpha = 0.2) +
    scale_color_manual(values = c(set_names(clone_colors, str_c(names(clone_colors), "_freq")), all_cells_freq = "black")) +
    scale_fill_manual(values = c(set_names(clone_colors, str_c(names(clone_colors), "_freq")), all_cells_freq = "black")) +
    scale_x_continuous(breaks = 1:15, labels = c(levels(cell_type_freq.tib$CellType), "")) +
    theme_bw() +
    xlab("") +
    ylab("Frequency") +
    theme(aspect.ratio = 1, panel.border = element_blank())

# Plot normalized frequency
pdf("210204_2_RadarPlot.pdf")
cell_type_freq2.tib %>% dplyr::select(CellType, contains(c("all_cells", current.clones))) %>%
    add_row(CellType = "empty") %>%
    replace(is.na(.), 0) %>%
    mutate(CellType = factor(CellType, levels = CellType)) %>%
    dplyr::select(CellType, contains("norm")) %>% pivot_longer(cols = -CellType) %>%
    mutate(name = factor(name, levels = str_c(current.clones, "_norm")[c(2:10,1)])) %>%
    ggplot(aes(x = as.numeric(CellType), y = value, fill = name, color = name)) +
    coord_radar() +
    geom_line() +
    geom_polygon(alpha = 0.2) +
    scale_color_manual(values = c(set_names(clone_colors, str_c(names(clone_colors), "_norm")), all_cells_norm = "black")) +
    scale_fill_manual(values = c(set_names(clone_colors, str_c(names(clone_colors), "_norm")), all_cells_norm = "black")) +
    scale_x_continuous(breaks = 1:15, labels = c(levels(cell_type_freq.tib$CellType), "")) +
    theme_bw() +
    xlab("") +
    ylab("Normalized cell type frequency") +
    theme(aspect.ratio = 1, panel.border = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "#D3D3D3"),
          axis.text = element_text(size = 14, color = "black"), axis.title.y = element_text(size = 14, color = "black"),
          axis.ticks.y = element_line(color = "black"))
dev.off()













