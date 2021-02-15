# Peter van Galen, 201101
# Make a UMAP of cell types in samples seu (diagnosis)
# This is based on 201101_SW_CellLineMix.R


#~~~~~~~~~~~~~~~~~~~~~~#
#### Pre-requisites ####
#~~~~~~~~~~~~~~~~~~~~~~#

options(stringsAsFactors = FALSE)
options(scipen = 999)

library(tidyverse)
library(data.table)
library(SingleCellExperiment)
library(Seurat)
library(gdata)

setwd("~/DropboxPartners/Projects/Maester/AnalysisPeter/210123_BPDCN712_Diagnosis")

rm(list=ls())

# Functions
source("../201007_FunctionsGeneral.R")
popcol.df <- read.xls("~/DropboxPartners/Pipelines/AuxiliaryFiles/PopCol.xlsx", sheet = 3, row.names = 1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Subset for common cells ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Load Seurat object from the Single-cell_BPDCN project
seu <- readRDS("~/DropboxPartners/Projects/Single-cell_BPDCN/AnalysisDaniel/201214_Seurat_GoT/BPDCN712_Seurat_Genotyped.rds")

# Load maegatk object
maegatk.rse <- readRDS("../200917_MT_Coverage/TenX_BPDCN712_mr3_maegatk.rds")

# Generate cell IDs that will match Maegtk
seu$cellMerge <- paste0(cutf(colnames(seu), d = "-", f = 1), "-1")

# Only keep cells with a cellMerge id that occurs once, then intersect with Maester data
cells.ch <- tibble(cell = seu$cellMerge) %>% group_by(cell) %>% filter(n()==1) %>% .$cell %>% unname
# 10095 / 10113 cells left (0.18% lost)
cells.ch <- cells.ch[cells.ch %in% colnames(maegatk.rse[[2]])]
# 9424 / 10095 cells left (6.6% lost)

# Subset for common cells
seu <- subset(seu, subset = cellMerge %in% cells.ch)
seu <- RenameCells(seu, new.names = as.character(seu$cellMerge))
seu$cellMerge <- NULL

# Exclude RPS and RPL genes
rp.ch <- rownames(seu)[grepl("^RPS|^RPL", rownames(seu))]
seu <- subset(seu, features = setdiff(rownames(seu), rp.ch))

# Normalize (log, transcript per 10K) and scale
seu <- NormalizeData(seu)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Add signatures of interest ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Define signatures of interest
BPDCN_Tumor.ch <- c("IGLL1", "PDLIM1", "TCL1A", "SERPINB2", "LAMP5", "HMSD", "PLVAP", "TNNI2", "STMN1", "SCN3A", "WASF1", "SERHL2", "NUCB2", "CLEC11A", "TPM2", "CEP70", "MYBPH", "SIRPG", "ARPP21", "LOXL4")
Cycle.ch <- na.omit(read.xls("~/DropboxPartners/Pipelines/AuxiliaryFiles/200727_signatures.xlsx")[-1,"tirosh_cycle"])

seu <- AddModuleScore(seu, features = list(BPDCN_Tumor.ch))
colnames(seu@meta.data) <- gsub("Cluster1", "BPDCN_Tumor", colnames(seu@meta.data))
seu <- AddModuleScore(seu, features = list(Cycle.ch))
colnames(seu@meta.data) <- gsub("Cluster1", "Cycle", colnames(seu@meta.data))
seu@meta.data %>% head


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Dimensionality reduction by marker genes ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# There is a single PreB cell that is actually a circulating tumor cells (see BPDCN paper); reclassify as pDC
plot(sort(seu$BPDCN_Tumor), col = ifelse(sort(seu$BPDCN_Tumor) > 0.75, "red", "black"), ylab = "Tumor score", xlab = "Cells")
subset(seu, subset = CellType == "PreB")@meta.data[,"BPDCN_Tumor"]
seu$CellType[seu$CellType == "PreB"] <- "pDC"
seu$CellType <- factor(seu$CellType, levels = levels(seu$CellType)[levels(seu$CellType) %in% unique(seu$CellType)])
seu@active.ident <- seu$CellType

# Define cluster marker genes (this takes a while, consider reading tsv at the end instead)
markerGenes <- FindAllMarkers(seu, slot = "data", logfc.threshold = 0.25, min.pct = 0.1, test.use = "roc", return.thresh = 0.4)
markergenes.dt.ls <- lapply(split(markerGenes, f = markerGenes$cluster), function(x) data.table(x))
markergenes.dt.ls <- lapply(markergenes.dt.ls, function(x) setorder(x, -avg_logFC))
# Check that all top 50 genes have a positive fold change
stopifnot(min(unlist(lapply(markergenes.dt.ls, function(x) x[1:50,avg_logFC]))) > 0)
markergenes.tib <- as_tibble( do.call(cbind, lapply(markergenes.dt.ls, function(x) x$gene[1:50])) )
write_tsv(markergenes.tib, file = "210123_MarkerGenes.txt")
markergenes.tib <- read_tsv("210123_MarkerGenes.txt")

# Scale data (don't use cell cycle regression which takes a long time and does not improve visualization)
#seu <- ScaleData(seu, features = rownames(seu), vars.to.regress = "Cycle")
seu <- ScaleData(seu)

# Linear dimension reduction. Use marker genes rather than variable genes to better match UMAP coordinates with clusters (that were determined with the random forest classification).
#seu <- RunPCA(seu, features = VariableFeatures(object = seu), seed.use = 42)
seu <- RunPCA(seu, features = unique(unname(unlist(markergenes.tib))), seed.use = 42)
DimPlot(seu, reduction = "pca", group.by = "replicate") +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    ggtitle("PCA colored by replicate")
DimHeatmap(seu, dims = 1:15, cells = 500, balanced = TRUE)

# Dimensionality reduction(s)
#for (ndims in c(5:20, 30, 40, 50)) {
ndims <- 9
seu <- RunUMAP(seu, dims = 1:ndims, seed.use = 42)

# Quick visualization
pdf(str_c("210123_1_UMAP_dim", ndims, ".pdf"))
print(
    DimPlot(seu, reduction = "umap", cols = popcol.df[levels(seu@active.ident),"hex"]) +
        theme(aspect.ratio = 1)
)
FeaturePlot(seu, features = "HBD")
FeaturePlot(seu, features = "Cycle")
dev.off()
#seu@reductions$umap <- NULL
#}

seu$UMAP_1 <- unname(seu@reductions$umap@cell.embeddings[,1])
seu$UMAP_2 <- unname(seu@reductions$umap@cell.embeddings[,2])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Exclude questionable cells ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# I also tried to run DecontX but the results do not make biological sense (e.g. all B cells 100% contaminated).

# Define markers genes for broad clusters (compare 201101_SW_CellLineMix.R)
seu$broad_cluster <- factor(gsub("Prog|ProMono|^Mono|ncMono|cDC|pDC", "Myeloid", gsub("HSC|EarlyE|LateE", "Erythroid", gsub("^B$|Plasma", "B", gsub("^T$|CTL|NK", "TNK", seu$CellType)))), levels = c("Myeloid", "Erythroid", "B", "TNK"))
seu$CellType %>% table
seu$broad_cluster %>% table
seu@active.ident <- seu$broad_cluster
broad_markers <- FindAllMarkers(seu, slot = "data", logfc.threshold = 0.25, min.pct = 0.1, test.use = "roc", return.thresh = 0.4)
broad_markers.dt.ls <- lapply(split(broad_markers, f = broad_markers$cluster), function(x) data.table(x))
broad_markers.dt.ls <- lapply(broad_markers.dt.ls, function(x) setorder(x, -avg_logFC))
# Check that all top 50 genes have a positive fold change
stopifnot(min(unlist(lapply(broad_markers.dt.ls, function(x) x[1:50,avg_logFC]))) > 0)
top_markers.tib <- as_tibble( do.call(cbind, lapply(broad_markers.dt.ls, function(x) x$gene[1:10])) )
seu@active.ident <- seu$CellType

# Calculate signature scores
for (n in colnames(top_markers.tib)) {
    seu <- AddModuleScore(object = seu, features = list(top_markers.tib[[n]]))
    colnames(seu@meta.data) <- gsub("Cluster1", n, colnames(seu@meta.data))
    }

# Visualize signature scores
for (n in colnames(top_markers.tib)) {
    print(
        FeaturePlot(seu, features = n, min.cutoff = "q25", order = T) +
            ggtitle(n) +
            theme(aspect.ratio = 1)
    )
    }

# Visualize scores for matching and non-matching broad cluster signatures
score_cutoff <- c(0, 0.5)
as_tibble(seu@meta.data, rownames = "cell") %>%
    pivot_longer(cols = colnames(top_markers.tib)) %>%
    ggplot() +
    geom_violin(aes(x = broad_cluster, y = value, fill = broad_cluster)) +
    geom_jitter(aes(x = broad_cluster, y = value), size = 0.2) +
    geom_hline(yintercept = score_cutoff, linetype = "dashed", col = "red") +
    facet_wrap(~name)

# Determine cells that have excessive expression of a marker genes that don't match their identity
contaminated.tib <- as_tibble(seu@meta.data, rownames = "cell") %>% select(cell, broad_cluster, Myeloid, Erythroid, B, TNK)
bad_cells_1 <- contaminated.tib %>% filter(B > 0.5, broad_cluster %in% c("Erythroid", "TNK"))
bad_cells_2 <- contaminated.tib %>% filter(Erythroid > 0.5, ! broad_cluster == "Erythroid")
bad_cells_3 <- contaminated.tib %>% filter(Myeloid > 0.5, broad_cluster == "TNK")
bad_cells_4 <- contaminated.tib %>% filter(TNK > 0, ! broad_cluster == "TNK")
bad_cells <- unique(c(bad_cells_1$cell, bad_cells_2$cell, bad_cells_3$cell, bad_cells_4$cell))
seu$contaminated <- colnames(seu) %in% bad_cells
#subset(seu, subset = contaminated == T)@meta.data %>% as_tibble(rownames = "cell") %>% view

pdf("210123_2_Contamination.pdf")
as_tibble(seu@meta.data, rownames = "cell") %>%
    pivot_longer(cols = colnames(top_markers.tib)) %>%
    ggplot() +
    geom_violin(aes(x = broad_cluster, y = value, fill = broad_cluster)) +
    geom_jitter(aes(x = broad_cluster, y = value, color = cell %in% bad_cells), size = 0.2) +
    scale_color_manual(values = c("black", "red")) +
    geom_hline(yintercept = score_cutoff, linetype = "dashed", col = "red") +
    facet_wrap(~name)

FeaturePlot(seu, features = "contaminated", cols = c("grey", "red"), order = T) +
    theme(aspect.ratio = 1)
dev.off()

# Exclude 78 / 9424 (0.8%) of cells based on contamination (leaving 9346)
seu <- subset(seu, contaminated == F)

# Plot
pdf(str_c("210123_3_UMAP_clean.pdf"), width = 6, height = 6)
print(
as_tibble(seu@meta.data, rownames = "cell") %>% #.$CellType %>% levels
    mutate(PlotOrder = factor(CellType, levels = c("CTL", "NK", "T", "EarlyE", "LateE",
                                                   "pDC", "cDC",
                                                   "ProMono", "Mono", "ncMono",
                                                   "Prog", "HSC",
                                                   "B", "Plasma"))) %>%
    arrange(PlotOrder) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = CellType)) +
    geom_point(size = 1) +
    scale_color_manual(values = popcol.df$hex) +
    theme_classic() +
    theme(aspect.ratio = 1, axis.line = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    guides(color = guide_legend(override.aes = list(size = 3)))
)
dev.off()


# Save results
as_tibble(seu@meta.data, rownames = "cell") %>% head
seu$broad_cluster <- NULL
seu$contaminated <- NULL
seu$orig.ident <- NULL
seu$tech <- NULL
seu@meta.data[,colnames(seu@meta.data) %in% colnames(top_markers.tib)] <- NULL
seu@meta.data[,grepl("Predict", colnames(seu@meta.data))] <- NULL
seu@meta.data[,grepl("tSNE", colnames(seu@meta.data))] <- NULL
seu <- DietSeurat(seu)
saveRDS(seu, file = "BPDCN712_Seurat.rds")

maegatk.rse <- maegatk.rse[[2]][,colnames(seu)]
saveRDS(maegatk.rse, file = "BPDCN712_Maegatk.rds")
