# Peter van Galen, 210219
# Process Seq-Well S^3 cell line mixing experiment data
# Goal is to save Seurat object with high-quality cells

#~~~~~~~~~~~~~~~~~~~~~~#
#### Pre-requisites ####
#~~~~~~~~~~~~~~~~~~~~~~#

options(stringsAsFactors = FALSE)
options(scipen = 999)

library(tidyverse)
library(readxl)
library(celda)
library(SingleCellExperiment)
library(Matrix)
library(Seurat)
library(data.table)
library(ggrastr)

setwd("~/DropboxPartners/Projects/Maester/AnalysisPeter/2_Cell_line_mixes_decontX")

rm(list=ls())

# Functions and colors (available at https://github.com/vangalenlab/MAESTER-2021)
source("../210215_FunctionsGeneral.R")
popcol.df <- read_excel("../MAESTER_colors.xlsx")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Cluster and marker genes ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Load expression data (available at https://vangalenlab.bwh.harvard.edu/maester-2021/)
load("../ExpressionData/SW_CellLineMix_expr.RData", verbose = T)
seu <- CreateSeuratObject(counts = CM.df, project = "SW_CellLineMix")
seu$Array <- cutf(colnames(CM.df), d = "\\.|_", f = 2)

# Generate cell IDs that will match MAEGATK output
seu$cellMerge <- paste0(cutf(colnames(seu), d = "_", f = 2), "-1")

# Only keep cells with a cellMerge id that occurs once
cell <- tibble(cell = seu$cellMerge) %>% group_by(cell) %>% filter(n()==1) %>% .$cell
seu <- seu[, seu$cellMerge %in% cell]

# Make the column names match MAEGATK & save
seu <- RenameCells(seu, new.names = as.character(seu$cellMerge))
seu$cellMerge <- NULL
write_tsv(tibble(sort(colnames(seu))), "SW_CellLineMix.QualityCells.txt", col_names = F)

# Exclude RPS and RPL genes
rp.ch <- rownames(seu)[grepl("^RPS|^RPL", rownames(seu))]
seu <- subset(seu, features = setdiff(rownames(seu), rp.ch))

# Normalize (log, transcript per 10K)
seu <- NormalizeData(seu)

# Identify variable features
seu <- FindVariableFeatures(seu)
LabelPoints(plot = VariableFeaturePlot(seu), points = head(VariableFeatures(seu), 10), repel = T, xnudge = 0, ynudge = 0) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + ggtitle("Variable genes")

# Scale data
seu <- ScaleData(seu, features = rownames(seu))

# Linear dimension reduction
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
DimPlot(seu, reduction = "pca", group.by = "Array") +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    ggtitle("PCA colored by array")
DimHeatmap(seu, dims = 1:9, cells = 500, balanced = TRUE)

# Cluster. Recommended resolution 0.4-1.2 for 3,000 cells, but we're doing broad clustering because we know there are two cell types.
seu <- FindNeighbors(seu, dims = 1:6)
seu <- FindClusters(seu, resolution = 0.05)

# Seurat visualizations
seu <- RunUMAP(seu, dims = 1:6)

# Quick visualizations to identify clusters
DimPlot(seu, reduction = "umap") + theme(aspect.ratio = 1)
FeaturePlot(seu, features = "HBG2", slot = "data", reduction = "umap") + theme(aspect.ratio = 1)
FeaturePlot(seu, features = "PTPRZ1", slot = "data", reduction = "umap") + theme(aspect.ratio = 1)
FeaturePlot(seu, features = "MKI67", slot = "data", reduction = "umap") + theme(aspect.ratio = 1)

# Give clusters intuitive names and save as CellType metadata
seu <- RenameIdents(object = seu, "0" = "K562", "1" = "BT142")
seu$CellType <- factor(seu@active.ident, levels = c("BT142", "K562"))
DimPlot(seu, reduction = "umap") + theme(aspect.ratio = 1)

# Also add UMAP coordinates to metadata
seu$UMAP_1 <- seu@reductions$umap@cell.embeddings[,1]
seu$UMAP_2 <- seu@reductions$umap@cell.embeddings[,2]

# Markers of gene expression using area under the curve
MarkerGenes.roc <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "roc")
write_tsv(MarkerGenes.roc, file = "SW_CellLineMix_ROC_Markers.txt")
MarkerGenes.roc <- read_tsv("SW_CellLineMix_ROC_Markers.txt")

# List of top 10 genes for each cluster
Top10Genes.ls <- lapply(split(MarkerGenes.roc, f = MarkerGenes.roc$cluster), function(x) x$gene[1:10])

# Add metadata and plot
pdf("SW_CellLineMix_1_Clusters_Scores.pdf", width = 6, height = 6)

DimPlot(seu, reduction = "umap") + theme(aspect.ratio = 1)

for (n in names(Top10Genes.ls)) {
    seu <- AddModuleScore(object = seu, features = Top10Genes.ls[n], name = n)
    colnames(seu@meta.data) <- gsub(str_c(n, "1$"), str_c(n, "_Score"), colnames(seu@meta.data))
    print(
        FeaturePlot(seu, features = paste0(n, "_Score"), min.cutoff = "q25", order = T) +
            ggtitle(n) +
            theme(aspect.ratio = 1)
    )
}

FeaturePlot(seu, features = "MKI67") + theme(aspect.ratio = 1)

dev.off()


#~~~~~~~~~~~~~~~#
#### DecontX ####
#~~~~~~~~~~~~~~~#

# Calculate contamination fraction & add to Seurat object
sce <- SingleCellExperiment(list(counts=as.matrix(GetAssayData(seu, slot = "counts"))))
sce <- decontX(x = sce, assayName = "counts") #, z = seu$CellType)
seu$decontX <- sce$decontX_contamination

# Visualize decontX results
pdf("SW_CellLineMix_2_DecontX.pdf", width = 6, height = 6)

# Color by contamination
cont_cutoff <- 0.05

# Color by contamination
ggplot(data = seu@meta.data, aes(x = UMAP_1, y = UMAP_2, fill = decontX, color = decontX > cont_cutoff)) + 
    geom_point(shape = 21) +
    scale_color_manual(values = c("NA", "red")) +
    scale_fill_gradientn(colors = c("blue","green","yellow","orange","red")) +
    theme_classic() +
    theme(aspect.ratio = 1)

# Histogram
ggplot(data = seu@meta.data, aes(x = decontX)) +
    geom_freqpoly(bins = 100) +
    geom_vline(xintercept = cont_cutoff, col = "red") +
    coord_cartesian(xlim=c(0,1)) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    ggtitle(str_c(round(mean(seu$decontX > cont_cutoff)*100,2), "% of cells have contamination > ", cont_cutoff))

dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Intersect module scores and decontX to exclude cells ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Plot evaluation of parameters to exclude cells
pdf("SW_CellLineMix_3_Cleanup.pdf", width = 6, height = 6)

# Cluster signature scores grouped by cluster
score_cutoff <- tibble(name = c("BT142_Score", "K562_Score"), cutoff = c(-0.5,0.2))
seu@meta.data %>% pivot_longer(cols = paste0(names(Top10Genes.ls), "_Score")) %>%
    ggplot(aes(x = CellType, y = value, fill = CellType)) +
    geom_violin() +
    geom_jitter(size = 0.2) +
    facet_wrap(~name) +
    geom_hline(data = score_cutoff, aes(yintercept = cutoff))

# Cells that are positive for either K562 or BT142
seu$MultiPositive <- apply(seu@meta.data[,c("BT142_Score", "K562_Score")], 1, function(x) sum(x > score_cutoff$cutoff) > 1)
table(seu$MultiPositive)

# Compare decontX to number of clusters
plot.dt <- data.table( seu@meta.data[,c("decontX", "MultiPositive")] )
setorder(plot.dt, -decontX)
par(mar=c(5,5,5,2))
plot(plot.dt$decontX, pch = 16, col = ifelse(plot.dt$MultiPositive, yes = NA, no = "black"), cex = 0.3,
     xlab = "cell", ylab = "deContX", main = "Red if scored high for >1 signature")
points(plot.dt$decontX, pch = 16, col = ifelse(plot.dt$MultiPositive, yes = "red", no = NA), cex = 0.3)
abline(h = cont_cutoff, col = "red")

# decontX score vs. number of clusters
VlnPlot(seu, features = "decontX", group.by = "MultiPositive", pt.size = 0.2) +
    labs(title = "Scores vs. decontX", x = "Positive score for >1 cluster", y = "Contamination estimate") +
    geom_hline(yintercept = cont_cutoff)

# Exclude if more than 1 *and* more than 10% contamination
seu$Keep <- ! (seu$MultiPositive & seu$decontX > cont_cutoff)
DimPlot(seu, group.by = "Keep", cols = c("grey", "red"), order = "FALSE") +
    ggtitle(str_c("Keep ", round(mean(seu$Keep)*100, 2), "% (exclude ",
                  round(100-mean(seu$Keep)*100,2), "%)")) +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))

dev.off()

# Exclude 218 / 2525 cells (8.63%) that did not pass; 2307 left
table(seu$Keep)
seu_red <- subset(seu, subset = Keep == T)
seu_red$CellType %>% table

# Save Seurat object
seu_red@commands <- list()
seu_red <- DietSeurat(seu_red, dimreducs = c("umap", "tsne"))
colnames(seu_red@meta.data)
seu_red$RNA_snn_res.0.05 <- NULL
seu_red$seurat_clusters <- NULL
seu_red$K562_Score <- NULL
seu_red$BT142_Score <- NULL
seu_red$MultiPositive <- NULL
seu_red$Keep <- NULL
saveRDS(seu_red, file = paste0(seu_red@project.name, "_Seurat_Keep.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Plot gene expression for Figure 1 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# You can skip everything above (except the prerequisites) and start here
heatcol.ch <- read_excel("../MAESTER_colors.xlsx", sheet = 2, col_names = "heatcol")$heatcol
seu_red <- readRDS(file = "SW_CellLineMix_Seurat_Keep.rds")
MarkerGenes.roc <- read_tsv("SW_CellLineMix_ROC_Markers.txt")
Top10Genes.ls <- lapply(split(MarkerGenes.roc, f = MarkerGenes.roc$cluster), function(x) x$gene[1:10])
genes.ch <- unlist(Top10Genes.ls)

expr.tib <- tibble(cell = colnames(seu_red), CellType = seu_red@meta.data$CellType,
                   UMAP_1 = seu_red$UMAP_1,
                   UMAP_2 = seu_red$UMAP_2,
                   data.frame( t(GetAssayData(seu_red, slot = "data")[genes.ch,])))

pdf("SW_CellLineMix_4_MarkerGenes.pdf")
for (g in genes.ch) {
print(
expr.tib %>%
    ggplot(aes_string(x = "UMAP_1", y = "UMAP_2", color = g)) +
    geom_point_rast() +
    scale_color_gradientn(colors = heatcol.ch[2:10]) +
    theme_classic() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    ggtitle(g)
)
}
dev.off()





