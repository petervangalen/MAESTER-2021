# Peter van Galen, 211012
# Trajectory analysis of erythroid and myeloid lineages

# Prerequisites
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(tidyverse)
library(Seurat)
library(slingshot)
library(SingleCellExperiment)
library(RColorBrewer)
library(readxl)
library(ggnewscale)
library(gridExtra)

rm(list=ls())
setwd("~/DropboxMGB/Projects/Maester/AnalysisPeter/7_Slingshot")

# Functions and colors (available at https://github.com/vangalenlab/MAESTER-2021)
source("../Auxiliary_files/210215_FunctionsGeneral.R")
popcol.df <- read_excel("../MAESTER_colors.xlsx")
mycol.ch <- popcol.df$hex
names(mycol.ch) <- popcol.df$name
timecol <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)

# Load Seurat object (available at https://vangalenlab.bwh.harvard.edu/maester-2021/))
seu <- readRDS("../4_CH_sample/BPDCN712_Seurat_Final.rds")
metadata.tib <- as_tibble(seu@meta.data, rownames = "cell")


# Trajectory --------------------------------------------------------------------------------------

# Function to get line x and y coordinates as well as pseudotime values for every cells
add_coordinates <- function(celltypes) {
    seu.subset <- subset(seu, subset = CellType %in% celltypes & UMAP_1 > -5) # UMAP filter is for a few distant pDCs
    
    # Slingshot
    sce <- as.SingleCellExperiment(seu.subset)
    reducedDim(sce, "umap") <- as.matrix(data.frame(seu.subset$UMAP_1, seu.subset$UMAP_2))
    sce <- slingshot(sce, clusterLabels = "CellType", reducedDim = "umap", start.clus = "HSC", approx_points = FALSE)
    
    # Visualize
    par(pty="s")
    plot(reducedDims(sce)$umap, col = timecol[cut(sce$slingPseudotime_1, breaks=100)], pch=16, asp = 1)
    lines(SlingshotDataSet(sce), lwd=2, col='black')
    
    # Extract trajectory
    sling_curve <- slingCurves(sce)[[1]]
    line_coord.df <- sling_curve[["s"]][sling_curve$ord,]
    line_coord.tib <- as_tibble(line_coord.df, rownames = "cell")
    line_coord.tib <- line_coord.tib %>%
        left_join(as_tibble(sce@colData[,"slingPseudotime_1",drop=F], rownames = "cell"), by = "cell")
    
    # Add trajectory info to Seurat metadata
    return(line_coord.tib)
    
}

### Erythroid: add trajectory and pseudotime to metadata
ery_coordinates.tib <- add_coordinates(celltypes = c("HSC", "EarlyE", "LateE"))
ery_coordinates.tib <- ery_coordinates.tib %>% rename(slingPseudotime_1 = "ery_pseudotime",
                                                      seu.subset.UMAP_1 = "ery_line_1", seu.subset.UMAP_2 = "ery_line_2")
metadata.tib <- metadata.tib %>% left_join(ery_coordinates.tib, by = "cell")

### Myeloid: add trajectory and pseudotime to metadata
mye_coordinates.tib <- add_coordinates(celltypes = c("HSC", "Prog", "ProMono", "Mono", "ncMono"))
mye_coordinates.tib <- mye_coordinates.tib %>% rename(slingPseudotime_1 = "mye_pseudotime",
                                                      seu.subset.UMAP_1 = "mye_line_1", seu.subset.UMAP_2 = "mye_line_2")
metadata.tib <- metadata.tib %>% left_join(mye_coordinates.tib, by = "cell")

### pDC: add trajectory and pseudotime to metadata
pdc_coordinates.tib <- add_coordinates(celltypes = c("HSC", "Prog", "pDC"))
pdc_coordinates.tib <- pdc_coordinates.tib %>% rename(slingPseudotime_1 = "pdc_pseudotime",
                                                      seu.subset.UMAP_1 = "pdc_line_1", seu.subset.UMAP_2 = "pdc_line_2")
metadata.tib <- metadata.tib %>% left_join(pdc_coordinates.tib, by = "cell")

### cDC: add trajectory and pseudotime to metadata
cdc_coordinates.tib <- add_coordinates(celltypes = c("HSC", "Prog", "cDC"))
cdc_coordinates.tib <- cdc_coordinates.tib %>% rename(slingPseudotime_1 = "cdc_pseudotime",
                                                      seu.subset.UMAP_1 = "cdc_line_1", seu.subset.UMAP_2 = "cdc_line_2")
metadata.tib <- metadata.tib %>% left_join(cdc_coordinates.tib, by = "cell")


### Visualize -------------------------------------------------------------------------------------

# Add clone info to metadata
positive_cells.tib <- read_tsv("../4_CH_sample/4.4_positive_cells.txt")
metadata.tib <- metadata.tib %>% mutate(clone = ifelse(cell %in% filter(positive_cells.tib, variant == "2593_G>A")$cell,
                                                       yes = "2593_G>A", no = "Other"))

# Plot each cell type on a separate page
celltype.ch <- c("mye", "ery", "cdc", "pdc")

pdf("7_Trajectories_Densities.pdf", width = 8, height = 5)

for (celltype in celltype.ch) {

# UMAP
g1 <- metadata.tib %>% arrange(get(str_c(celltype, "_pseudotime"))) %>%
    filter(CellType %in% c("HSC", "Prog", "EarlyE", "LateE", "ProMono", "Mono", "ncMono", "cDC", "pDC") & UMAP_1 > -5) %>%
    ggplot(mapping = aes_string(x = str_c(celltype, "_line_1"), y = str_c(celltype, "_line_2"))) +
    geom_point(mapping = aes_string(x = "UMAP_1", y = "UMAP_2", fill = str_c(celltype, "_pseudotime"),
                                    color = "clone"), size = 2, pch = 21, stroke = 0.25) +
    scale_fill_gradientn(colors = timecol, na.value = "#DDDDDD") +
    scale_color_manual(values = c("black", "#FFFFFF00")) +
    geom_path(col = "black", size = 0.3) +
    ggtitle(celltype) +
    theme_classic() +
    theme(aspect.ratio = 1,
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())

# Density
g2 <- metadata.tib %>%
    ggplot(aes_string(x = str_c(celltype, "_pseudotime"))) +
    geom_density(aes(color = clone), adjust = 1.5) +
    scale_color_manual(values = c("#7b68ee", "#48d1cc")) +
    new_scale_color() +
    geom_point(data = metadata.tib, aes(y = ifelse(clone == "Other", 0, 0.03), color = CellType), size = 0.5) +
    scale_color_manual(values = mycol.ch[levels(metadata.tib$CellType)]) +
    ylab("density") +
    theme_classic()+
    theme(aspect.ratio = 0.75,
          axis.text = element_text(color = "black"),
          axis.ticks = element_line(color = "black"))

grid.arrange(g1, g2, layout_matrix = matrix(c(1,2), nrow = 1))

}

dev.off()


