# Peter van Galen, 210204
# Assess lineage bias in the clonal hematopoiesis sample, including the radar plot in Supplemental Figure 16


#~~~~~~~~~~~~~~~~~~~~~#
#### Prerequisites ####
#~~~~~~~~~~~~~~~~~~~~~#

options(scipen = 999)

library(tidyverse)
library(Seurat)
library(SummarizedExperiment)
library(Matrix)
library(ggrastr)
library(ComplexHeatmap)
library(readxl)

rm(list=ls())
setwd("~/DropboxMGB/Projects/Maester/AnalysisPeter/4_CH_sample")

# Functions and colors (available on https://github.com/vangalenlab/MAESTER-2021)
source("../210215_FunctionsGeneral.R")
popcol.df <- read_excel("../MAESTER_colors.xlsx")
mycol.ch <- popcol.df$hex
names(mycol.ch) <- popcol.df$name

# Load Seurat object (available at https://vangalenlab.bwh.harvard.edu/maester-2021/)
seu <- readRDS("BPDCN712_Seurat_Final.rds")
metadata.tib <- as_tibble(seu@meta.data, rownames = "cell")

# Load Maegtk (also (available at https://vangalenlab.bwh.harvard.edu/maester-2021/), then calculate allele frequencies
maegatk.rse <- readRDS("BPDCN712_Maegatk_Final.rds")
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))*100
# Check (should be TRUE)
all(colnames(af.dm) == colnames(seu))

# Variants of interest
voi.ch <- read.table("4.3_vois_order.txt")$V1


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

# Save a tibble with cell, variant, coverage, clone, cell type information
positive_metadata.tib <- as_tibble(bind_rows(lapply(positive_cells.ls, function(x) data.frame(cell = x)), .id = "variant")[,2:1]) %>%
    mutate(variant = factor(variant, levels = voi.ch))
cov.mat <- as.matrix( assays(maegatk.rse)$coverage )
positive_metadata.tib$cov <- apply(positive_metadata.tib, 1, function(x) { cov.mat[as.numeric(cutf(x[2], d = "_")),x[1]] } )
# Add clone information (see also "Group similar variants" in "4.3_Variants_Of_Interest.R")
positive_metadata.tib <- positive_metadata.tib %>%
    mutate(clone = gsub("6205_G>A|9164_T>C", "6205_G>A-9164_T>C", gsub(
        "10158_T>A|6293_T>C", "10158_T>A-6293_T>C", gsub(
            "1415_G>A|9753_G>A", "1415_G>A-9753_G>A", variant)))) %>%
    mutate(clone = factor(clone, levels = unique(clone)))
# Add cell type information
positive_metadata.tib <- positive_metadata.tib %>% left_join(dplyr::select(as_tibble(seu@meta.data, rownames = "cell"), cell, CellType))
write_tsv(positive_metadata.tib, file = "4.4_positive_cells.txt")

# Save a tibble with cell type counts per clone
cell_type_numbers.tib <- metadata.tib %>% group_by(CellType) %>% dplyr::count(name = "all_cells") %>% ungroup
cell_type_numbers_variants.tib <- positive_metadata.tib %>% select(variant, CellType) %>% group_by(variant, CellType) %>% summarize(number = n()) %>%
    pivot_wider(id_cols = CellType, names_from = variant, values_from = number) %>% replace(is.na(.), 0)
cell_type_numbers.tib <- left_join(cell_type_numbers.tib, cell_type_numbers_variants.tib, by = "CellType")
write_tsv(cell_type_numbers.tib, file = "4.4_cell_type_numbers.txt")


#~~~~~~~~~~~~~~~~#
#### Bar plot ####
#~~~~~~~~~~~~~~~~#

# Transpose, pivot
transpose_numbers.tib <- cell_type_numbers.tib %>% pivot_longer(c(all_cells, all_of(voi.ch)), names_to = "clone") %>% pivot_wider(names_from = CellType)

plot.tib <- transpose_numbers.tib %>%
    mutate(clone = factor(clone, levels = c("all_cells", setdiff(clone, "all_cells")))) %>%
    pivot_longer(-clone, names_to = "CellType", values_to = "Number of cells") %>%
    mutate(CellType = factor(CellType, levels = popcol.df$name[1:14]))

pdf("4.4_1_CellTypeBarplots.pdf")

plot.tib %>% filter(clone != "all_cells") %>%
ggplot(aes(x = clone, y = `Number of cells`, fill = CellType)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = mycol.ch[levels(plot.tib$CellType)]) +
    scale_y_continuous(limits=c(0, 300)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 0.7)

ggplot(plot.tib, aes(x = clone, y = `Number of cells`, fill = CellType)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = mycol.ch[levels(plot.tib$CellType)]) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 0.7)

dev.off()

# Add frequency and normalized frequency
cell_type_freq.tib <- cell_type_numbers.tib %>%
    mutate(across(-CellType, .fns = list(freq = ~ .x / sum(.x) * 100, # equivalent to function(x) {x / sum(x) * 100 }
                                         norm = ~ (.x / sum(.x)) / (all_cells / sum(all_cells)))))

# Save for Tyler's P-value calcuations
cell_type_freq.tib %>% dplyr::select(CellType, contains("freq")) %>% write_tsv(file = "4.4_cell_type_freq.txt")
cell_type_freq.tib %>% dplyr::select(CellType, contains("norm")) %>% write_tsv(file = "4.4_cell_type_norm.txt")


#~~~~~~~~~~~~~~~~~~#
#### Radar plot ####
#~~~~~~~~~~~~~~~~~~#

# Now, use clones instead of variants
cell_type_numbers2.tib <- metadata.tib %>% group_by(CellType) %>% dplyr::count(name = "all_cells") %>% ungroup
cell_type_numbers_clones.tib <- positive_metadata.tib %>% filter(!duplicated(cell)) %>%
    select(clone, CellType) %>% group_by(clone, CellType) %>% summarize(number = n()) %>%
    pivot_wider(id_cols = CellType, names_from = clone, values_from = number) %>% replace(is.na(.), 0)
cell_type_numbers2.tib <- left_join(cell_type_numbers2.tib, cell_type_numbers_clones.tib, by = "CellType")

# Add frequency and normalized frequency
cell_type_freq2.tib <- cell_type_numbers2.tib %>%
    mutate(across(-CellType, .fns = list(freq = ~ .x / sum(.x) * 100, # equivalent to function(x) {x / sum(x) * 100 }
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
    scale_color_manual(values = c(set_names(mycol.ch, str_c(names(mycol.ch), "_freq")), all_cells_freq = "black")) +
    scale_fill_manual(values = c(set_names(mycol.ch, str_c(names(mycol.ch), "_freq")), all_cells_freq = "black")) +
    scale_x_continuous(breaks = 1:15, labels = c(levels(cell_type_freq.tib$CellType), "")) +
    theme_bw() +
    xlab("") +
    ylab("Frequency") +
    theme(aspect.ratio = 1, panel.border = element_blank())

# Plot normalized frequency
radar.tib <- cell_type_freq2.tib %>% dplyr::select(CellType, contains(c("all_cells", current.clones))) %>%
    add_row(CellType = "empty") %>%
    replace(is.na(.), 0) %>%
    mutate(CellType = factor(CellType, levels = CellType)) %>%
    dplyr::select(CellType, contains("norm")) %>% pivot_longer(cols = -CellType) %>%
    mutate(name = factor(name, levels = str_c(current.clones, "_norm")[c(2:10,1)]))
# Prepare colors
clone_colors <- mycol.ch[intersect(names(mycol.ch), gsub("_norm", "", levels(radar.tib$name)))]
clone_colors <- c(set_names(clone_colors, str_c(names(clone_colors), "_norm")), all_cells_norm = "black")

pdf("4.4_2_RadarPlot.pdf")
radar.tib  %>%
    ggplot(aes(x = as.numeric(CellType), y = value, fill = name, color = name)) +
    coord_radar() +
    geom_line() +
    geom_polygon(alpha = 0.2) +
    scale_color_manual(values = clone_colors) +
    scale_fill_manual(values = clone_colors) +
    scale_x_continuous(breaks = 1:15, labels = c(levels(cell_type_freq.tib$CellType), "")) +
    theme_bw() +
    xlab("") +
    ylab("Normalized cell type frequency") +
    theme(aspect.ratio = 1,
          panel.border = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "#D3D3D3"),
          axis.text = element_text(size = 14, color = "black"), axis.title.y = element_text(size = 14, color = "black"),
          axis.ticks.y = element_line(color = "black"))
dev.off()

# In the final figure, only control, 2593G>A and 6243G>A are shown (others are removed in Illustrator)











