# Peter van Galen, 210123
# Import single-cell analysis of an individual with clonal hematopoesis; make a UMAP save a new Seurat object
# This sample is also referred to as Patient 10 (Griffin et al., under review) and BPDCN712


#~~~~~~~~~~~~~~~~~~~~~~#
#### Pre-requisites ####
#~~~~~~~~~~~~~~~~~~~~~~#

options(stringsAsFactors = FALSE)
options(scipen = 999)

library(tidyverse)
library(data.table)
library(SingleCellExperiment)
library(Seurat)
library(readxl)

setwd("~/DropboxPartners/Projects/Maester/AnalysisPeter/4_CH_sample")

rm(list=ls())

# Functions
source("../210215_FunctionsGeneral.R")
popcol.df <- read_excel("../MAESTER_colors.xlsx")
mycol.ch <- popcol.df$hex
names(mycol.ch) <- popcol.df$name


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Subset for common cells ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Load Seurat object from the clonal hematopoiesis bone marrow aspirate (available at https://vangalenlab.bwh.harvard.edu/resources/single-cell_bpdcn/)
seu <- readRDS("~/DropboxPartners/Projects/Single-cell_BPDCN/AnalysisDaniel/201214_Seurat_GoT/BPDCN712_Seurat_Genotyped.rds")

# Load maegatk object (available at https://vangalenlab.bwh.harvard.edu/maester-2021/)
se.ls <- readRDS("../1_MT_Coverage/TenX_BPDCN712_mr3_maegatk.rds")
maegatk.rse <- se.ls[[2]]

# Generate cell IDs that will match Maegtk
seu$cellMerge <- paste0(cutf(colnames(seu), d = "-", f = 1), "-1")

# Only keep cells with a cellMerge id that occurs once, then intersect with Maester data
cells.ch <- tibble(cell = seu$cellMerge) %>% group_by(cell) %>% filter(n()==1) %>% .$cell %>% unname
# 10095 / 10113 cells left (0.18% lost)
cells.ch <- cells.ch[cells.ch %in% colnames(maegatk.rse)]
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

# Define signatures of interest: tumor cell signature from Griffin et al., under review
BPDCN_Tumor_score.ch <- c("IGLL1", "PDLIM1", "TCL1A", "SERPINB2", "LAMP5", "HMSD", "PLVAP", "TNNI2", "STMN1", "SCN3A", "WASF1", "SERHL2", "NUCB2", "CLEC11A", "TPM2", "CEP70", "MYBPH", "SIRPG", "ARPP21", "LOXL4")
# And from Tirosh et al., https://www.nature.com/articles/nature20123, Table S1
Cycle.ch <- c("HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA", "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8")

seu <- AddModuleScore(seu, features = list(BPDCN_Tumor_score.ch))
colnames(seu@meta.data) <- gsub("Cluster1", "BPDCN_Tumor_score", colnames(seu@meta.data))
seu <- AddModuleScore(seu, features = list(Cycle.ch))
colnames(seu@meta.data) <- gsub("Cluster1", "Cycle_score", colnames(seu@meta.data))
seu@meta.data %>% head


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Dimensionality reduction by marker genes ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# There is a single PreB cell that is actually a circulating tumor cells (see BPDCN paper); reclassify as pDC
seu$CellType %>% table
subset(seu, subset = CellType == "PreB")@meta.data[,"BPDCN_Tumor_score"]
plot(sort(seu$BPDCN_Tumor_score), col = ifelse(sort(seu$BPDCN_Tumor_score) > 0.75, "red", "black"), ylab = "Tumor score", xlab = "Cells")
seu$CellType[seu$CellType == "PreB"] <- "pDC"
seu$CellType <- factor(seu$CellType, levels = levels(seu$CellType)[levels(seu$CellType) %in% seu$CellType])
seu@active.ident <- seu$CellType

# Define cluster marker genes (this takes a while, consider reading tsv at the end instead)
markerGenes <- FindAllMarkers(seu, slot = "data", logfc.threshold = 0.25, min.pct = 0.1, test.use = "roc", return.thresh = 0.4, only.pos = T)
markergenes.dt.ls <- lapply(split(markerGenes, f = markerGenes$cluster), function(x) data.table(x))
markergenes.dt.ls <- lapply(markergenes.dt.ls, function(x) setorder(x, -avg_logFC))
markergenes.tib <- as_tibble( do.call(cbind, lapply(markergenes.dt.ls, function(x) x$gene[1:50])) )
write_tsv(markergenes.tib, file = "210123_MarkerGenes.txt")
markergenes.tib <- read_tsv("210123_MarkerGenes.txt")

# Scale data (don't use cell cycle regression which takes a long time and does not improve visualization)
#seu <- ScaleData(seu, features = rownames(seu), vars.to.regress = "Cycle")
seu <- ScaleData(seu)

# Linear dimension reduction. Use marker genes rather than variable genes to better match UMAP coordinates with clusters (that were determined with the random forest classification; Griffin et al., under review).
#seu <- RunPCA(seu, features = VariableFeatures(object = seu), seed.use = 42)
seu <- RunPCA(seu, features = unique(unname(unlist(markergenes.tib))), seed.use = 42)
DimPlot(seu, reduction = "pca", group.by = "replicate") +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    ggtitle("PCA colored by replicate")
DimHeatmap(seu, dims = 1:15, cells = 500, balanced = TRUE)

# Dimensionality reduction. I settled on 9 dimensions empirically.
ndims <- 9
seu <- RunUMAP(seu, dims = 1:ndims, seed.use = 42)

# Quick visualization
pdf(str_c("210123_1_UMAP_dim", ndims, ".pdf"))
print(
    DimPlot(seu, reduction = "umap", cols = mycol.ch) +
        theme(aspect.ratio = 1)
)
FeaturePlot(seu, features = "HBD")
FeaturePlot(seu, features = "Cycle_score")
dev.off()

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
broad_markers <- FindAllMarkers(seu, slot = "data", logfc.threshold = 0.25, min.pct = 0.1, test.use = "roc", return.thresh = 0.4, only.pos = T)
broad_markers.dt.ls <- lapply(split(broad_markers, f = broad_markers$cluster), function(x) data.table(x))
broad_markers.dt.ls <- lapply(broad_markers.dt.ls, function(x) setorder(x, -avg_logFC))
top_markers.tib <- as_tibble( do.call(cbind, lapply(broad_markers.dt.ls, function(x) x$gene[1:10])) )
seu@active.ident <- seu$CellType

# Calculate signature scores
for (n in colnames(top_markers.tib)) {
    seu <- AddModuleScore(object = seu, features = list(top_markers.tib[[n]]))
    colnames(seu@meta.data) <- gsub("Cluster1", str_c(n, "_score"), colnames(seu@meta.data))
    }

# Visualize signature scores
for (n in str_c(colnames(top_markers.tib), "_score")) {
    print(
        FeaturePlot(seu, features = n, min.cutoff = "q25", order = T) +
            ggtitle(n) +
            theme(aspect.ratio = 1)
    )
    }

# Visualize scores for matching and non-matching broad cluster signatures
score_cutoff <- c(0, 0.5)
as_tibble(seu@meta.data, rownames = "cell") %>%
    pivot_longer(cols = str_c(colnames(top_markers.tib), "_score")) %>%
    ggplot() +
    geom_violin(aes(x = broad_cluster, y = value, fill = broad_cluster)) +
    geom_jitter(aes(x = broad_cluster, y = value), size = 0.2) +
    geom_hline(yintercept = score_cutoff, linetype = "dashed", col = "red") +
    facet_wrap(~name)

# Determine cells that have excessive expression of a marker genes that don't match their identity
contaminated.tib <- as_tibble(seu@meta.data, rownames = "cell") %>% select(cell, CellType, broad_cluster, Myeloid_score,
                                                                           Erythroid_score, B_score, TNK_score)
bad_cells_1 <- contaminated.tib %>% filter(B_score > 0.5, broad_cluster %in% c("Erythroid", "TNK"))
bad_cells_2 <- contaminated.tib %>% filter(Erythroid_score > 0.5, ! broad_cluster == "Erythroid")
bad_cells_3 <- contaminated.tib %>% filter(Myeloid_score > 0.5, broad_cluster == "TNK")
bad_cells_4 <- contaminated.tib %>% filter(TNK_score > 0, ! broad_cluster == "TNK")
bad_cells <- unique(c(bad_cells_1$cell, bad_cells_2$cell, bad_cells_3$cell, bad_cells_4$cell))
seu$contaminated <- colnames(seu) %in% bad_cells
#subset(seu, subset = contaminated == T)@meta.data %>% as_tibble(rownames = "cell") %>% view

pdf("210123_2_Contamination.pdf")
as_tibble(seu@meta.data, rownames = "cell") %>%
    pivot_longer(cols = str_c(colnames(top_markers.tib), "_score")) %>%
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
    scale_color_manual(values = mycol.ch) +
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
seu@meta.data[,colnames(seu@meta.data) %in% str_c(colnames(top_markers.tib), "_score")] <- NULL
seu@meta.data[,grepl("Predict", colnames(seu@meta.data))] <- NULL
seu@meta.data[,grepl("tSNE", colnames(seu@meta.data))] <- NULL
seu <- DietSeurat(seu)
saveRDS(seu, file = "BPDCN712_Seurat.rds")

maegatk.rse <- maegatk.rse[,colnames(seu)]
saveRDS(maegatk.rse, file = "BPDCN712_Maegatk.rds")
