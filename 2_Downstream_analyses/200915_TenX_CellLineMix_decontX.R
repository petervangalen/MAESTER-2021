# Peter van Galen, 200915
# Process 10X V3 Cell line mixing experiment data
# Goal is to save Seurat object with cells that have low contamination

options(max.print = 500)
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(tidyverse)
library(data.table)
library(celda)
library(SingleCellExperiment)
library(Matrix)
library(Seurat)
library(gdata)
library(ggplot2)

setwd("~/DropboxPartners/Projects/Maester/AnalysisPeter/200915_All_Clustering_Decontx")

rm(list=ls())

# Functions
source("../201007_FunctionsGeneral.R")
popcol.df <- read.xls("~/DropboxPartners/Pipelines/AuxiliaryFiles/PopCol.xlsx", sheet = 3, row.names = 1)


#===========================================
# Broad clusters and marker genes
#===========================================

# This section resembles "200819_AML_clustering.R" and the first part of "200820_AML_cleanup.R"
load("../CellRanger/K562-BT142_agg/K562-BT142_agg.RData", verbose = T)
seu <- CreateSeuratObject(counts = CM.dgm, project = "TenX_CellLineMix")
seu$Lane <- cutf(colnames(seu), d = "-", f = 2)

# Generate cell IDs that will match Maegtk
seu$cellMerge <- paste0(cutf(colnames(seu), d = "-"), "-1")

# Only keep cells with a cellMerge id that occurs once
cell <- tibble(cell = seu$cellMerge) %>% group_by(cell) %>% filter(n()==1) %>% .$cell
seu <- seu[, seu$cellMerge %in% cell]

# Make the column names match Maegtk
seu <- RenameCells(seu, new.names = as.character(seu$cellMerge))
seu$cellMerge <- NULL

# Save QualityCells for Maegtk filter
write_tsv(data.frame(colnames(seu)), "TenX_CellLineMix.QualityCells.txt", col_names = F)

# Normalize (log, transcript per 10K)
seu <- NormalizeData(seu)

# Identify variable features
seu <- FindVariableFeatures(seu)
LabelPoints(plot = VariableFeaturePlot(seu), points = head(VariableFeatures(seu), 10), repel = T, xnudge = 0, ynudge = 0) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + ggtitle("Variable genes")

# Scale data
seu <- ScaleData(seu, features = rownames(seu))

# Linear dimension reduction
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
DimPlot(seu, reduction = "pca", group.by = "Lane") +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    ggtitle("PCA colored by lane")
DimHeatmap(seu, dims = 1:9, cells = 500, balanced = TRUE)

# Cluster. Recommended resolution 0.4-1.2 for 3,000 cells, but we're doing broad clustering because we know there are two cell types.
seu <- FindNeighbors(seu, dims = 1:6)
seu <- FindClusters(seu, resolution = 0.1)

# Seurat visualizations
seu <- RunUMAP(seu, dims = 1:6)
seu <- RunTSNE(seu, dims = 1:6)

# Quick visualizations to identify clusters
pdf("TenX_CellLineMix_1_Clusters.pdf", width = 6, height = 6)
DimPlot(seu, reduction = "umap") + theme(aspect.ratio = 1)
DimPlot(seu, reduction = "tsne") + theme(aspect.ratio = 1)
FeaturePlot(seu, features = "HBG2", slot = "data", reduction = "umap") + theme(aspect.ratio = 1)
FeaturePlot(seu, features = "PTPRZ1", slot = "data", reduction = "umap") + theme(aspect.ratio = 1)
FeaturePlot(seu, features = "MKI67", slot = "data", reduction = "umap") + theme(aspect.ratio = 1)
dev.off()

# Give clusters intuitive names and save as CellType metadata
seu <- RenameIdents(object = seu, "0" = "K562", "1" = "BT142", "2" = "BT142_Cycling")
seu$CellType <- factor(seu@active.ident, levels = c("BT142", "K562", "BT142_Cycling"))
DimPlot(seu, reduction = "umap") + theme(aspect.ratio = 1)

# Markers of gene expression using area under the curve
MarkerGenes.roc <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "roc")
write.table(MarkerGenes.roc, file = "TenX_CellLineMix_ROC_Markers.txt", sep = "\t" , quote = FALSE, col.names = NA)


#===========================================
# DecontX
#===========================================

# Calculate contamination fraction & add to Seurat object
sce <- SingleCellExperiment(list(counts=as.matrix(GetAssayData(seu, slot = "counts"))))
sce <- decontX(x = sce, assayName = "counts") #, z = seu$CellType)
seu$decontX <- sce$decontX_contamination

# Quick visualizations
pdf("TenX_CellLineMix_2_DecontX.pdf", width = 6, height = 6)
umap <- reducedDim(sce, "decontX_UMAP")
plotDimReduceCluster(umap[,1], umap[,2], cluster = sce$decontX_clusters)
plotDecontXContamination(sce)
plotDecontXContamination(sce, colorScale = c(rep("green",100), rep("red",900)))
ggplot(data = seu@meta.data, aes(x = decontX)) +
    geom_freqpoly(bins = 100) +
    geom_vline(xintercept = 0.1, col = "red") +
    coord_cartesian(xlim=c(0,1))
dev.off()


#===========================================
# Add module scores
#===========================================

# List of top 10 genes for each cluster, determined by ROC above
Top10Genes.ls <- lapply(split(MarkerGenes.roc, f = MarkerGenes.roc$cluster), function(x) x$gene[1:10])

# Add metadata and plot
pdf("TenX_CellLineMix_3_Scores.pdf", width = 6, height = 6)
for (n in names(Top10Genes.ls)) {
    seu <- AddModuleScore(object = seu, features = Top10Genes.ls[n], name = paste0(n, "_"))
    print(
        FeaturePlot(seu, features = paste0(n, "_1"), min.cutoff = "q25", order = T) +
            ggtitle(n) +
            theme(aspect.ratio = 1)
    )
}
dev.off()


# Plot evaluation of parameters to exclude cells
pdf("TenX_CellLineMix_4_Contamination.pdf", width = 6, height = 6)

# Plot populations to orient
DimPlot(seu, group.by = "CellType", pt.size = 0.3) +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + ggtitle("Named clusters")

# Cluster signature scores grouped by cluster
VlnPlot(seu, features = paste0(names(Top10Genes.ls), "_1"), pt.size = 0.3)

# Contamination score
FeaturePlot(seu, features = "decontX", cols = c("blue", "green", "yellow", "orange", "red"))

# Cells that are positive for either K562 or BT142
seu$SinglePositive <- apply(seu@meta.data[,c("K562_1", "BT142_1")], 1, function(x) sum(x > 0) == 1)
table(seu$SinglePositive)

# Compare decontX to number of clusters
plot.dt <- data.table( seu@meta.data[,c("decontX", "SinglePositive")] )
setorder(plot.dt, -decontX)
par(mar=c(5,5,5,2))
plot(plot.dt$decontX, pch = 16, col = ifelse(plot.dt$SinglePositive, yes = "black", no = "NA"), cex = 0.3,
     xlab = "cell", ylab = "deContX", main = "Red if scored higher than 0 for both cell types")
points(plot.dt$decontX, pch = 16, col = ifelse(plot.dt$SinglePositive, yes = NA, no = "red"), cex = 0.3)
abline(h = 0.1, col = "red")

# decontX score vs. number of clusters
VlnPlot(seu, features = "decontX", group.by = "SinglePositive", pt.size = 0.1) +
    labs(title = "Scores vs. decontX", x = "Positive score for exactly one cluster", y = "Contamination estimate")

# Exclude if more than 1 *and* more than 10% contamination
seu$Keep <- ! (!seu$SinglePositive & seu$decontX > 0.1)
DimPlot(seu, group.by = "Keep", cols = c("grey", "red"), order = "FALSE") + ggtitle("Keep") +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))

dev.off()

# How many cells are you excluding with these parameters?
table(seu$Keep)

# Save Seurat object
seu@commands <- list()
seu <- DietSeurat(seu, dimreducs = c("umap", "tsne"))
colnames(seu@meta.data)
seu$RNA_snn_res.0.1 <- NULL
seu$seurat_clusters <- NULL
seu$K562_1 <- NULL
seu$BT142_1 <- NULL
seu$BT142_Cycling_1 <- NULL
seu$SinglePositive <- NULL
saveRDS(seu, file = paste0(seu@project.name, "_Seurat_Keep.rds"))





