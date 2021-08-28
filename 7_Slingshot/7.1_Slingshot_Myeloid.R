# Peter van Galen, 210209
# Trajectory analysis of myeloid lineage

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

rm(list=ls())
setwd("~/DropboxPartners/Projects/Maester/AnalysisPeter/7_Slingshot")

# Functions and colors (available on https://github.com/vangalenlab/MAESTER-2021)
source("../210215_FunctionsGeneral.R")
popcol.df <- read_excel("../MAESTER_colors.xlsx")
mycol.ch <- popcol.df$hex
names(mycol.ch) <- popcol.df$name

# Load Seurat object (available at https://vangalenlab.bwh.harvard.edu/maester-2021/)
seu <- readRDS("../4_CH_sample/BPDCN712_Seurat_Final.rds")


#~~~~~~~~~~~~~~~~~~#
#### Trajectory ####
#~~~~~~~~~~~~~~~~~~#

# What cells to select?
seu.subset1 <- subset(seu, subset = CellType %in% c("HSC", "Prog", "ProMono", "Mono", "ncMono"))

# Slingshot
#rm(sling)
sling <- as.SingleCellExperiment(seu.subset1)
reducedDim(sling, "umap") <- as.matrix(data.frame(seu.subset1$UMAP_1, seu.subset1$UMAP_2))
sling <- slingshot(sling, clusterLabels = "CellType", reducedDim = "umap")

# Visualize
par(pty="s")
timecol <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plot(reducedDims(sling)$umap, col = timecol[cut(sling$slingPseudotime_1, breaks=100)], pch=16, asp = 1)
lines(SlingshotDataSet(sling), lwd=2, col='black')

# Extract trajectory
SlingCurve <- slingCurves(sling)[[1]]
SlingCurve <- SlingCurve[["s"]][SlingCurve$ord,]
SlingCurve <- tibble(cell = rownames(SlingCurve),
                     SlingOrder = 1:nrow(SlingCurve),
                     SlingCurve_1 = SlingCurve[,1],
                     SlingCurve_2 = SlingCurve[,2])
SlingCurve <- SlingCurve %>% left_join(as_tibble(sling@colData[,"slingPseudotime_1",drop=F], rownames = "cell"), by = "cell") %>%
    rename(slingPseudotime_1 = "Pseudotime")

# Add trajectory info to Seurat metadata
metadata.tib <- as_tibble(seu@meta.data, rownames = "cell")
metadata.tib <- metadata.tib %>% left_join(SlingCurve, by = "cell")

# Add clone info to metadata
positive_cells.tib <- read_tsv("../4_CH_sample/4.3_positive_cells.txt")
metadata.tib <- metadata.tib %>% mutate(clone = ifelse(cell %in% filter(positive_cells.tib, variant == "2593_G>A")$cell,
                                                       yes = "2593_G>A", no = "Other"))

pdf("7.1_1_Trajectory_analysis.pdf")
print(
metadata.tib %>% arrange(SlingOrder) %>%
    filter(CellType %in% c("HSC", "Prog", "EarlyE", "LateE", "ProMono", "Mono", "ncMono", "cDC", "pDC") & UMAP_1 > -5) %>%
    ggplot(mapping = aes(x = SlingCurve_1, y = SlingCurve_2)) +
    geom_point(mapping = aes(x = UMAP_1, y = UMAP_2, fill = Pseudotime, color = clone), size = 3, pch = 21) +
    scale_fill_gradientn(colors = timecol, na.value = "#DDDDDD") +
    scale_color_manual(values = c("black", "#FFFFFF00")) +
    geom_path(col = "black", size = 1) +
    theme_classic() +
    theme(aspect.ratio = 1, axis.line = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
)
dev.off()

# Add pseudotime to Seurat meta.data & save
stopifnot(all(metadata.tib$cell == colnames(seu)))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Density over pseudotime ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Plot density across trajectory
pdf("7.1_2_Density.pdf")
metadata.tib %>%
    ggplot(aes(x = Pseudotime)) +
    geom_density(aes(color = clone)) +
    scale_color_manual(values = c("#7b68ee", "#48d1cc")) +
    new_scale_color() +
    geom_point(data = metadata.tib, aes(y = ifelse(clone == "Other", 0, 0.03), color = CellType)) +
    scale_color_manual(values = mycol.ch) +
    scale_x_continuous(breaks = c(0, 2.5, 5)) +
    theme_classic()+
    theme(aspect.ratio = 0.75)
dev.off()





