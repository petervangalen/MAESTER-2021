# Peter van Galen, 201119
# Detect subclones in K562 cells from the SW cell line mixing experiment
# This is based on 201118_SW_K562_clones.R


#~~~~~~~~~~~~~~~~~~~~~#
#### Prerequisites ####
#~~~~~~~~~~~~~~~~~~~~~#

options(stringsAsFactors = FALSE)
options(scipen = 999)

library(tidyverse)
library(magrittr)
library(SummarizedExperiment)
library(Seurat)
library(data.table)
library(Matrix)
library(ComplexHeatmap)
library(gdata)
library(ggrastr)
library(ggforce)
#library(stringr)

setwd("~/DropboxPartners/Projects/Maester/AnalysisPeter/200922_Cell_line_mixes")
rm(list=ls())

# More functions
source("../201007_FunctionsGeneral.R")
popcol.df <- read.xls("~/DropboxPartners/Pipelines/AuxiliaryFiles/PopCol.xlsx", sheet = 3, row.names = 1)
heatcol.ch <- read.xls("~/DropboxPartners/Pipelines/AuxiliaryFiles/PopCol.xlsx", sheet = 4, header = F)$V1

### Import data for Seq-Well S^3 Cell Line Mix
# Maegatk output
se.ls <- readRDS(paste0("../200917_MT_Coverage/SW_CellLineMix_mr3_maegatk.rds"))
# Don't need Seurat object:
#seu <- readRDS("../200915_All_Clustering_Decontx/SW_CellLineMix_Seurat_Keep.rds")
# This tibble contains UMAP coordinates and cell type classifications (see 201101_SW_CellLineMix.R)
cells.tib <- read_rds("SW_CellLineMix_Cells.rds")

# For 10X 3' v3 Cell Line Mix
#seu <- readRDS("../200915_All_Clustering_Decontx/TenX_CellLineMix_Seurat_Keep.rds")
#se.ls <- readRDS(paste0("../200917_MT_Coverage/TenX_CellLineMix_mr3_Maegtk.rds"))


#~~~~~~~~~~~~~~~~~~~~#
#### Prepare data ####
#~~~~~~~~~~~~~~~~~~~~#

# Put maegatk object in the same order as cells.tib
maegatk.rse <- se.ls[[2]][,cells.tib$cell]

# Prepare allele frequency matrix
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))*100

# Group cell IDs by cell type classification. Use the cell type classification from maegatk, not RNA-seq.
table( as.character(cells.tib$CellType_RNA) == as.character(cells.tib$CellType_MT) )
CellSubsets.ls <- list(unionCells = cells.tib$cell,
                       K562 = filter(cells.tib, CellType_MT == "K562")$cell,
                       BT142 = filter(cells.tib, CellType_MT == "BT142")$cell)

# Check coverage between cell types
pdf("SW_Coverage_between_CellTypes.pdf")
cells.tib %>% mutate(coverage = maegatk.rse$depth) %>%
    ggplot(aes(x = CellType_RNA, y = coverage)) +
    geom_violin() +
    theme(aspect.ratio = 1)
cells.tib %>% mutate(coverage = maegatk.rse$depth) %>% arrange(-coverage) %>% mutate(key = row_number()) %>%
    ggplot(aes(x = key, y = coverage, color = CellType_RNA)) +
    geom_point() +
    xlab("Rank sorted cells")
mean( sort(maegatk.rse$depth, decreasing = T)[1:500] )
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Collect variant information ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# This section takes a while; consider skipping to read_tsv at the end.

# Get the mean allele frequency and coverage for every cell subset.
mean_af.ls <- lapply(CellSubsets.ls, function(x) rowMeans(af.dm[,x]))
mean_cov.ls <- lapply(CellSubsets.ls, function(x) rowMeans(assays(se.ls[[2]])[["coverage"]][,x])[as.numeric(cutf(rownames(af.dm), d = "_"))])
names(mean_af.ls) <- paste0("mean_af.", names(mean_af.ls))
names(mean_cov.ls) <- paste0("mean_cov.", names(mean_cov.ls))

# Get the quantiles of the VAFs of each variant in each cell subset. This takes 2.5 minutes.
quantiles <- c("q01" = 0.01, "q05" = 0.05, "q10" = 0.1, "q50" = 0.5, "q90" = 0.9, "q95" = 0.95, "q99" = 0.99)
# This can take a while.
start_time <- Sys.time()
quantiles.ls <- lapply(quantiles, function(x) lapply(CellSubsets.ls, function(y) apply(af.dm[,y], 1, quantile, x) ))
Sys.time() - start_time

# Get the mean quality for each variant. This takes an hour.
start_time <- Sys.time()
qual.num <- sapply(rownames(af.dm), function(x) {
    pos <- as.numeric( cutf(x, d = "_") )
    message(pos)
    mut <- cutf(x, d = ">", f = 2)
    # Only use cells in which the base was sequenced. Use reverse only because that's how we amplify transcripts.
    covered <- assays(se.ls[[2]])[[str_c(mut, "_counts_rev")]][pos,] > 0
    # Get the mean quality for this call
    qual <- mean( assays(se.ls[[2]])[[str_c(mut, "_qual_rev")]][pos, covered] )
    return(qual)
})
Sys.time() - start_time

# Collect all information in a tibble
var.tib <- as_tibble(do.call(cbind, c(mean_af.ls, mean_cov.ls, unlist(quantiles.ls, recursive = F))), rownames = "var")
var.tib <- add_column(var.tib, quality = qual.num, .before = 2)

# Plot quality metrics
var.tib %>% arrange(quality) %>% mutate(key = row_number()) %>%
    mutate(mut = cutf(var, d = ">", f = 2)) %>%
    ggplot(aes(x = key, y = quality, color = mut)) +
    geom_point() +
    xlab(label = "Rank sorted variant") +
    theme(aspect.ratio = 0.5)

# Save for fast loading next time
write_tsv(var.tib, "SW_CellLineMix_Variants2.txt")
var.tib <- read_tsv("SW_CellLineMix_Variants2.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Detect variants of interest ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

### Get the most interesting subclonal variants
voi.ch <- var.tib %>%
    filter(mean_cov.unionCells > 20, q01.K562 < 1, q99.K562 > 2,
           mean_af.K562 > (5*mean_af.BT142),
           quality >= 30) %>%
    .$var
# Replace 9117_T>C, which comes up because it has no coverage in some cells, with another homoplasmic variant (2141_T>C), that has better coverage (see also 201101_SW_CellLineMix.R)
voi.ch <- str_replace_all(voi.ch, "9117_T>C", "2141_T>C")
str_view(voi.ch, "G>A|A>G|C>T|T>C")
#filter(var.tib, var %in% voi.ch) %>% view

# Order from large to small clones.
clone_sizes <- sapply(voi.ch, function(x) sum(af.dm[x,] > 1) )
voi.ch <- names(clone_sizes[rev(order(clone_sizes))])

# Add coverage and allele frequency info from variants of interest to cells.tib
cells.tib <- read_rds("SW_CellLineMix_Cells.rds") # same as above
for (voi in voi.ch) {
    pos <- as.numeric(cutf(voi, d = "_"))
    wt <- cutf(voi, d = "_|>|\\.", f = 2)
    mut <- cutf(voi, d = ">|\\.", f = 2)
    cells.tib <- cells.tib %>%
        left_join(as_tibble(assays(maegatk.rse)[["coverage"]][pos,], rownames = "cell"), by = "cell") %>%
        left_join(as_tibble(assays(maegatk.rse)[[str_c(mut, "_qual_rev")]][pos,], rownames = "cell"), by = "cell") %>%
        left_join(as_tibble(af.dm[voi,], rownames = "cell"), by = "cell") %>%
        rename(value.x = str_c("cov_", str_replace(voi, ">", ".")),
               value.y = str_c("qual_rev_", str_replace(voi, ">", ".")),
               value = str_c("af_", str_replace(voi, ">", ".")))
}
#cells.tib %>% filter(CellType_MT == "K562") %>% ggplot(aes(x = af_9117_T.C, y = cov_9117_T.C)) + geom_point()

### Plot VAF in each cell, sorted from low to high, illustrating the selection process
pdf("SW_K562_clones_1_Sorted_VAFs.pdf", width = 10, height = 5)
par(mar=c(4,4,2,8), xpd = T)
mycol <- rainbow(length(voi.ch))
plot(NA, xlim = c(0, nrow(cells.tib)), ylim = c(0, 100), xlab = "Cells (sorted separately for each variant)", ylab = "Variant Allele Frequency")
for (n in 1:length(voi.ch)) {
    v <- str_c("af_", str_replace(voi.ch[n], ">", "."))
    points(as.numeric(sort(cells.tib[[v]])), pch = 16, col = mycol[n])
    text(x = 1.1*nrow(cells.tib), y = 110-4*n, label = v, col = mycol[n], pos = 4, cex = 0.8)
}
dev.off()

### For each variant of interest, plot UMAP of cells, colored by VAF.
pdf(paste0("SW_K562_clones_2_VAF_UMAP.pdf"))
for (voi in voi.ch) {
    cov_colname <- str_c("cov_", str_replace(voi, ">", "."))
    af_colname <- str_c("af_", str_replace(voi, ">", "."))
    
    # Plot VAF on UMAP coordinates
    print(
        cells.tib %>% # filter(.[[cov_colname]] > 3) %>%  # optional: select cells with three genotyped transcripts
            arrange(.[[af_colname]]) %>%
            ggplot(aes_string(x = "UMAP_1", y = "UMAP_2", color = af_colname)) +
            geom_point_rast() +
            scale_color_gradientn(colors = heatcol.ch[2:10]) +
            theme_classic() + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
            ggtitle(voi)
    )
}
dev.off()


#~~~~~~~~~~~~~~~~#
#### Heatmaps ####
#~~~~~~~~~~~~~~~~#

### Prepare matrix of variants of interest in K562 cells, capped at VAF 15%
af_subset.mat <- af.dm[voi.ch,filter(cells.tib, CellType_MT == "K562")$cell]
af_subset.mat[af_subset.mat > 15] <- 15

### Variant correlation. Exclude the first two that are very different
var_cor.mat <- af_subset.mat[! rownames(af_subset.mat) %in% c("2141_T>C", "8213_G>A"),]
cor.mat <- cor(t(var_cor.mat))
# Cluster
var.clust <- hclust(as.dist(1-cor.mat))
ngroups <- 15

# Correlation heatmap
hm1 <- Heatmap(cor.mat,
               cluster_columns = var.clust,
               cluster_rows = var.clust,
               row_split = ngroups,
               column_split = ngroups,
               row_gap = unit(1, "mm"),
               column_gap = unit(1, "mm"),
               width = unit(100, "mm"),
               height = unit(100, "mm"),
               column_title = ngroups)
pdf(paste0("SW_K562_clones_3_Correlation.pdf"))
print(hm1)
dev.off()

# Determine variants in each group as displayed by Heatmap package:
Groups.tib <- tibble(var = names(cutree(var.clust, k = ngroups)), Cut = cutree(var.clust, k = ngroups))[var.clust$order,]
Groups.tib <- Groups.tib %>% mutate(Group = match(Cut, unique(Cut)))
Groups.ls <- split(Groups.tib$var, Groups.tib$Group)
# Order from Groups with most cells to groups least cells
GroupSizes.ls <- lapply(Groups.ls, function(x) sum(af.dm[x,] > 1))
Groups.ls <- Groups.ls[rev(order(unlist(GroupSizes.ls)))]
# And then from Groups with many to few variants
Groups.ls <- Groups.ls[order(-lengths(Groups.ls))]

### Customize column order
# First, sort for the large variant 8213_G>A
plot_order.mat <- af_subset.mat[,order(-af_subset.mat["8213_G>A",])]
# Then, sort for all variants from the correlation matrix
plot_order.mat <- plot_order.mat[unlist(Groups.ls),]
# Only use VAFs of >2% for sorting
plot_order.mat[plot_order.mat < 2 ] <- 0
# Order from high to low
for (x in rev(lapply(Groups.ls, function(x) match(x, unlist(Groups.ls))))) {
    if (length(x) == 1) { # for individual variants, order by VAF
        plot_order.mat <- plot_order.mat[,order(-plot_order.mat[x,])]
    } else { # for variants in a Group, use Pearson clustering for ordering
        tmp.mat <- plot_order.mat[x,! colSums(plot_order.mat[x,]) == 0]
        tmp.cor <- cor(tmp.mat, method = "pearson")
        tmp.hclust <- hclust(as.dist(1-tmp.cor))
        tmp_order.mat <- plot_order.mat[,tmp.hclust$labels[tmp.hclust$order]]
        plot_order.mat <- cbind(tmp_order.mat, plot_order.mat[,!colnames(plot_order.mat) %in% colnames(tmp_order.mat)])
    }
}

# Generate a matrix to plot, while maintaining sorting order, adding back large variants, and restoring VAFs of 0-2%.
plot.mat <- af_subset.mat[c("2141_T>C", "8213_G>A", rownames(plot_order.mat)), colnames(plot_order.mat)]

# Remove uninformative cells
plot_reduced.mat <- plot.mat[,colSums(plot.mat[! rownames(plot.mat) %in% c("2141_T>C", "8213_G>A"),] > 2) | plot.mat["8213_G>A",] < 15]
dim(plot.mat); dim(plot_reduced.mat)

# Plot
hm2 <- Heatmap(plot_reduced.mat,
               col = c("#FCFCFC","#FFEDB0","#FFDF5F","#FEC510","#FA8E24","#F14C2B","#DA2828","#BE2222","#A31D1D"),
               show_row_names = T,
               show_column_names = F,
               cluster_columns = F,
               cluster_rows = F,
               name = "AF",
               heatmap_legend_param = list(at = c(0, 5, 10, 15), border = "#000000", grid_height = unit(10, "mm")),
               border = T,
               width = unit(100, "mm"),
               height = unit(100, "mm"),
               use_raster = T)
pdf(paste0("SW_K562_clones_4_Heatmap.pdf"))
print(hm2)
dev.off()


















