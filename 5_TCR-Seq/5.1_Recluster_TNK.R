# Peter van Galen, 210827
# Recluster T/NK cells to get more granularity

# Prerequisites
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(tidyverse)
library(Seurat)
library(readxl)
library(harmony)
library(gridExtra)
library(data.table)
library(viridis)

setwd("~/DropboxMGB/Projects/Maester/AnalysisPeter/5_TCR-Seq/")

rm(list=ls())

# Functions and colors (available at https://github.com/vangalenlab/MAESTER-2021)
source("../210215_FunctionsGeneral.R")
popcol.tib <- read_excel("../MAESTER_colors.xlsx")
mycol.ch <- popcol.tib$hex
names(mycol.ch) <- popcol.tib$name

# Load Seurat object (available at https://vangalenlab.bwh.harvard.edu/maester-2021/)
seu <- readRDS("../4_CH_sample/BPDCN712_Seurat_Final.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Dimensionality reduction and clustering ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Select a subset of cells
tnk <- subset(seu, CellType %in% c("T", "CTL", "NK"))

# Scale data
tnk <- ScaleData(tnk)

# Linear dimension reduction by initial PCA
tnk <- FindVariableFeatures(tnk)
tnk <- RunPCA(tnk, features = VariableFeatures(tnk), seed.use = 42)
DimPlot(tnk, reduction = "pca", group.by = "replicate") +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    ggtitle("Initial PCA")

# Integrate to remove batch effect
tnk <- RunHarmony(object = tnk, group.by.vars = "replicate", reduction = "pca", plot_convergence = T)
DimPlot(tnk, reduction = "harmony", group.by = "replicate") + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    ggtitle("Harmony PCA")

# Cluster & visualize results
pdf("5.1.1_ReclusterTNK.pdf", width = 8, height = 8)

# Dimensionality reduction & clustering. I settled on 9 dimensions empirically.
#all_dims <- 5:10
#for (ndims in all_dims) {
ndims <- 9
message(ndims)
tnk <- RunUMAP(tnk, reduction = "harmony", dims = 1:ndims, seed.use = 42)
tnk <- FindNeighbors(tnk, reduction = "harmony", dims = 1:ndims)
tnk <- FindClusters(tnk, resolution = 0.4)

# Generate informative plots
rpl <- DimPlot(tnk, reduction = "umap", group.by = "replicate") +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    ggtitle(str_c("Replicate, ", ndims, " dim"))

celltype <- DimPlot(tnk, reduction = "umap", group.by = "CellType", cols = mycol.ch) +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    ggtitle("Group by cell type")

cycleScore <- FeaturePlot(tnk, features = "Cycle_score") +
    scale_color_viridis() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    ggtitle("Cell cycle score")

clust <- DimPlot(tnk, reduction = "umap") +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    ggtitle(str_c("Clusters, ", ndims, " dim"))

grid.arrange(rpl, celltype, cycleScore, clust, layout_matrix = matrix(c(1:4), nrow = 2, byrow = T))

#}
dev.off()

# Define cluster marker genes (this takes a while, consider reading tsv at the end instead)
markerGenes <- FindAllMarkers(tnk, slot = "data", logfc.threshold = 0.25, min.pct = 0.1,
                              test.use = "roc", return.thresh = 0.4, only.pos = T)
markergenes.dt.ls <- lapply(split(markerGenes, f = markerGenes$cluster), function(x) data.table(x))
markergenes.dt.ls <- lapply(markergenes.dt.ls, function(x) setorder(x, -avg_logFC))
markergenes.tib <- as_tibble( do.call(cbind, lapply(markergenes.dt.ls, function(x) x$gene[1:50])) )

# Save UMAP coordinates
tnk$UMAP_TNK_1 <- unname(tnk@reductions$umap@cell.embeddings[,1])
tnk$UMAP_TNK_2 <- unname(tnk@reductions$umap@cell.embeddings[,2])

# Annotate clusters based on Tyler and Zeyu's annotations
tnk <- RenameIdents(object = tnk,
                    "0" = "Hyperactivated Effector T cells",
                    "1" = "Effector T cells",
                    "2" = "Central memory T cells 1",
                    "3" = "Central memory T cells 2",
                    "4" = "Effector Memory T cells",
                    "5" = "NK cells",
                    "6" = "Early activated T cells",
                    "7" = "Proliferating lymphocyte")

# Save the marker genes
colnames(markergenes.tib) <- c("Hyperactivated Effector T cells", "Effector T cells", "Central memory T cells 1",
                               "Central memory T cells 2", "Effector Memory T cells", "NK cells",
                               "Early activated T cells", "Proliferating lymphocyte")
markergenes.tib <- markergenes.tib[,c(3,4,7,5,2,1,6,8)]
write_tsv(markergenes.tib, file = "TNK_MarkerGenes.txt")

# Wrangle
tnk <- DietSeurat(tnk)
tnk@meta.data$RNA_snn_res.0.4 <- NULL
tnk@meta.data$seurat_clusters <- NULL
tnk@meta.data$TNK_CellType <- factor(tnk@active.ident, levels = c("Central memory T cells 1", "Central memory T cells 2",
    "Early activated T cells", "Effector Memory T cells", "Effector T cells", "Hyperactivated Effector T cells",
    "NK cells", "Proliferating lymphocyte"))

pdf("5.1.2_AnnotatedTNK.pdf")

as_tibble(tnk@meta.data, rownames = "cell") %>%
    ggplot(aes(x = UMAP_TNK_1, y = UMAP_TNK_2, color = TNK_CellType)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = mycol.ch) +
    theme_bw() +
    theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 3) ) )

dev.off()

saveRDS(tnk, file = "TNK_Seurat.rds")
