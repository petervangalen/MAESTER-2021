# Peter van Galen, 201101
# Process Seq-Well S^3 Cell line mixing experiment data
# Goal is to save Seurat object with cells that have low contamination

#~~~~~~~~~~~~~~~~~~~~~~#
#### Pre-requisites ####
#~~~~~~~~~~~~~~~~~~~~~~#

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
library(RColorBrewer)
library(ggforce)
library(ggrastr)

setwd("~/DropboxPartners/Projects/Maester/AnalysisPeter/200915_All_Clustering_Decontx")

rm(list=ls())

# Functions
source("../201007_FunctionsGeneral.R")
popcol.df <- read.xls("~/DropboxPartners/Pipelines/AuxiliaryFiles/PopCol.xlsx", sheet = 3, row.names = 1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Cluster and marker genes ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# This section resembles "200819_AML_clustering.R" and the first part of "200820_AML_cleanup.R"
load("../ExprStar/190909.190830.K562BT142.Array12/2K.UMI-1K.gene-0.2rRNA-0.2chrM/190909.190830.K562BT142.Array12.filter.RData", verbose = T)
seu <- CreateSeuratObject(counts = CM.df, project = "SW_CellLineMix")
seu$Array <- cutf(colnames(CM.df), d = "\\.|_", f = 2)

# Generate cell IDs that will match Maegtk
seu$cellMerge <- paste0(cutf(colnames(seu), d = "_", f = 2), "-1")

# Only keep cells with a cellMerge id that occurs once
cell <- tibble(cell = seu$cellMerge) %>% group_by(cell) %>% filter(n()==1) %>% .$cell
seu <- seu[, seu$cellMerge %in% cell]

# Make the column names match Maegtk & save
seu <- RenameCells(seu, new.names = as.character(seu$cellMerge))
seu$cellMerge <- NULL
write_tsv(tibble(sort(colnames(seu))), "SW_CellLineMix.QualityCells.txt", col_names = F)

# Exclude RPS and RPL genes that are different between the cell lines but often obscure known markers.
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
seu <- RunTSNE(seu, dims = 1:6)

# Quick visualizations to identify clusters
pdf("SW_CellLineMix_1_Clusters.pdf", width = 6, height = 6)
DimPlot(seu, reduction = "umap") + theme(aspect.ratio = 1)
DimPlot(seu, reduction = "tsne") + theme(aspect.ratio = 1)
FeaturePlot(seu, features = "HBG2", slot = "data", reduction = "umap") + theme(aspect.ratio = 1)
FeaturePlot(seu, features = "PTPRZ1", slot = "data", reduction = "umap") + theme(aspect.ratio = 1)
FeaturePlot(seu, features = "MKI67", slot = "data", reduction = "umap") + theme(aspect.ratio = 1)
dev.off()

# Give clusters intuitive names and save as CellType metadata
seu <- RenameIdents(object = seu, "0" = "K562", "1" = "BT142")
seu$CellType <- factor(seu@active.ident, levels = c("BT142", "K562"))
DimPlot(seu, reduction = "umap") + theme(aspect.ratio = 1)

# Markers of gene expression using area under the curve
MarkerGenes.roc <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.score_cutoff = 0.25, test.use = "roc")
write.table(MarkerGenes.roc, file = "SW_CellLineMix_ROC_Markers.txt", sep = "\t" , quote = FALSE, col.names = NA)


#~~~~~~~~~~~~~~~#
#### DecontX ####
#~~~~~~~~~~~~~~~#

# Calculate contamination fraction & add to Seurat object
sce <- SingleCellExperiment(list(counts=as.matrix(GetAssayData(seu, slot = "counts"))))
sce <- decontX(x = sce, assayName = "counts") #, z = seu$CellType)
seu$decontX <- sce$decontX_contamination

# Quick visualizations
pdf("SW_CellLineMix_2_DecontX.pdf", width = 6, height = 6)

# Color by contamination
cont_cutoff <- 0.1
umap <- reducedDim(sce, "decontX_UMAP")
plotDimReduceCluster(umap[,1], umap[,2], cluster = sce$decontX_clusters)
plotDecontXContamination(sce)
plotDecontXContamination(sce, colorScale = c(rep("green",cont_cutoff*1000), rep("red",(1-cont_cutoff)*1000)))

# Histogram
ggplot(data = seu@meta.data, aes(x = decontX)) +
    geom_freqpoly(bins = 100) +
    geom_vline(xintercept = 0.1, col = "red") +
    coord_cartesian(xlim=c(0,1)) +
    theme_classic() +
    theme(aspect.ratio = 1)

# Color by contamination differently
CM.stats$CellMerge <- paste0(cutf(CM.stats$cell, d = "_", f = 2), "-1")
CM1.stats <- CM.stats[match(colnames(seu), CM.stats$CellMerge)]
ggplot(CM1.stats, aes(x = tSNEx, y = tSNEy, fill = seu$decontX, color = seu$decontX > 0.1)) + 
    geom_point(shape = 21) +
    scale_color_manual(values = c("NA", "red")) +
    theme_classic() +
    theme(aspect.ratio = 1)

dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Add module scores ####
#~~~~~~~~~~~~~~~~~~~~~~~~~#

# List of top 10 genes for each cluster, determined by ROC above
Top10Genes.ls <- lapply(split(MarkerGenes.roc, f = MarkerGenes.roc$cluster), function(x) x$gene[1:10])

# Add metadata and plot
pdf("SW_CellLineMix_3_Scores.pdf", width = 6, height = 6)
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
pdf("SW_CellLineMix_4_Contamination.pdf", width = 6, height = 6)

# Plot populations to orient
DimPlot(seu, group.by = "CellType", pt.size = 0.3) +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + ggtitle("Named clusters")

# Contamination score
mycol <- rev( colorRampPalette(brewer.pal(11,'Spectral')[-6])(100) )
tibble(cbind(seu@meta.data, seu@reductions$umap@cell.embeddings)) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, fill = decontX, color = seu$decontX > cont_cutoff)) +
    geom_point(shape = 21) +
    scale_color_manual(values = c("NA", "red")) +
    scale_fill_gradientn(colors = mycol) +
    theme_classic() +
    theme(aspect.ratio = 1)

# Cluster signature scores grouped by cluster
score_cutoff <- 0.2
tibble(cbind(seu@meta.data, seu@reductions$umap@cell.embeddings)) %>%
    pivot_longer(cols = paste0(names(Top10Genes.ls), "_1")) %>%
    ggplot() +
    geom_violin(aes(x = CellType, y = value, fill = CellType)) +
    geom_jitter(aes(x = CellType, y = value), size = 0.2) +
    geom_hline(yintercept = score_cutoff) +
    facet_wrap(~name)

# Cells that are positive for either K562 or BT142
seu$MultiPositive <- apply(seu@meta.data[,c("K562_1", "BT142_1")], 1, function(x) sum(x > score_cutoff) > 1)
table(seu$MultiPositive)

# Compare decontX to number of clusters
plot.dt <- data.table( seu@meta.data[,c("decontX", "MultiPositive")] )
setorder(plot.dt, -decontX)
par(mar=c(5,5,5,2))
plot(plot.dt$decontX, pch = 16, col = ifelse(plot.dt$MultiPositive, yes = NA, no = "black"), cex = 0.3,
     xlab = "cell", ylab = "deContX", main = paste0("Red if scored higher than ", score_cutoff, " for >1 clustertypes"))
points(plot.dt$decontX, pch = 16, col = ifelse(plot.dt$MultiPositive, yes = "red", no = NA), cex = 0.3)
abline(h = cont_cutoff, col = "red")

# decontX score vs. number of clusters
VlnPlot(seu, features = "decontX", group.by = "MultiPositive", pt.size = 0.2) +
    labs(title = "Scores vs. decontX", x = "Positive score for >1 cluster", y = "Contamination estimate") +
    geom_hline(yintercept = cont_cutoff)

# Exclude if more than 1 *and* more than 10% contamination
seu$Keep <- ! (seu$MultiPositive & seu$decontX > cont_cutoff)
DimPlot(seu, group.by = "Keep", cols = c("grey", "red"), order = "FALSE") + ggtitle("Keep") +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))

dev.off()

# Exclude 140 / 2525 cells (5.54%) that did not pass; 2385 left
table(seu$Keep)
seu_red <- subset(seu, subset = Keep == T)

# Save Seurat object
seu_red@commands <- list()
seu_red <- DietSeurat(seu_red, dimreducs = c("umap", "tsne"))
colnames(seu_red@meta.data)
seu_red$RNA_snn_res.0.05 <- NULL
seu_red$seurat_clusters <- NULL
seu_red$K562_1 <- NULL
seu_red$BT142_1 <- NULL
seu_red$MultiPositive <- NULL
seu_red$Keep <- NULL
saveRDS(seu_red, file = paste0(seu_red@project.name, "_Seurat_Keep.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Plot gene expression for Figure 1 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# You can skip everything above (except the prerequisites) and start here
rm(list=ls())
heatcol.ch <- read.xls("~/DropboxPartners/Pipelines/AuxiliaryFiles/PopCol.xlsx", sheet = 4, header = F)$V1
seu_red <- readRDS(file = "SW_CellLineMix_Seurat_Keep.rds")
MarkerGenes.roc <- read.table(file = "SW_CellLineMix_ROC_Markers.txt")
Top10Genes.ls <- lapply(split(MarkerGenes.roc, f = MarkerGenes.roc$cluster), function(x) x$gene[1:10])
genes.ch <- unlist(Top10Genes.ls)

expr.tib <- tibble(cell = colnames(seu_red), CellType = seu_red@meta.data$CellType,
                   UMAP_1 = seu_red@reductions$umap@cell.embeddings[,"UMAP_1"],
                   UMAP_2 = seu_red@reductions$umap@cell.embeddings[,"UMAP_2"],
                   data.frame( t(GetAssayData(seu_red, slot = "data")[genes.ch,])))

pdf("SW_CellLineMix_5_MarkerGenes.pdf")
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





