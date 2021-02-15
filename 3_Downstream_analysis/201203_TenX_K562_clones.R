# Peter van Galen, 201203
# Detect subclones in K562 cells from the 10X cell line mixing experiment
# This is based on 201119_SW_K562_clones.R


#~~~~~~~~~~~~~~~~~~~~~#
#### Prerequisites ####
#~~~~~~~~~~~~~~~~~~~~~#

options(stringsAsFactors = FALSE)
options(scipen = 999)

library(tidyverse)
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

    
#~~~~~~~~~~~~~~~~~~~~#
#### Prepare data ####
#~~~~~~~~~~~~~~~~~~~~#

# This tibble contains UMAP coordinates and cell type classifications (see 201101_SW_CellLineMix.R)
cells.tib <- read_rds("TenX_CellLineMix_Cells.rds")

# This tibble contains information about variants
vars.tib <- read_tsv("TenX_CellLineMix_Variants1.txt")

# Additional data can be accessed through the Maegatk output and Seurat object:
se.ls <- readRDS(paste0("../200917_MT_Coverage/TenX_CellLineMix_Combined_mr3_maegatk.rds"))
maegatk.rse <- se.ls[[2]][,cells.tib$cell]
#seu <- readRDS("../200915_All_Clustering_Decontx/TenX_CellLineMix_Seurat_Keep.rds")

# Prepare allele frequency matrix
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))*100

# Group cell IDs by cell type classification. Use the cell type classification from maegatk, not RNA-seq.
cells.tib <- cells.tib %>% mutate(CellType_RNA = gsub("BT142_Cycling", "BT142", CellType_RNA))
cells.tib$CellType_RNA %>% table
cells.tib$CellType_MT %>% table
table( as.character(cells.tib$CellType_RNA) == as.character(cells.tib$CellType_MT) )

# Check coverage between cell types
pdf("TenX_Coverage_between_CellTypes.pdf")
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
#### Detect variants of interest ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

### Get the most interesting subclonal variants
voi.ch <- vars.tib %>%
    filter(mean_cov.K562 > 100, q01.K562 < 1, q99.K562 > 2,
           mean_af.K562 > (5*mean_af.BT142),
           quality >= 30) %>%
    .$var
# Replace 9117_T>C, which comes up because it has no coverage in some cells, with another homoplasmic variant (2141_T>C), that has better coverage (see also 201101_SW_CellLineMix.R)
#voi.ch <- str_replace_all(voi.ch, "9117_T>C", "2141_T>C")
#str_view(voi.ch, "G>A|A>G|C>T|T>C")
#filter(vars.tib, var %in% voi.ch) %>% view

# Order from large to small clones.
clone_sizes <- sapply(voi.ch, function(x) sum(af.dm[x,] > 1) )
voi.ch <- names(clone_sizes[rev(order(clone_sizes))])

# Add coverage and allele frequency info from variants of interest to cells.tib
cells.tib <- read_rds("TenX_CellLineMix_Cells.rds") # same as above
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

### Plot VAF in each cell, sorted from low to high, illustrating the selection process
pdf("TenX_K562_clones_1_Sorted_VAFs.pdf", width = 10, height = 5)
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
pdf(paste0("TenX_K562_clones_2_VAF_UMAP.pdf"))
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

### Variant correlation & cluster
cor.mat <- cor(t(af_subset.mat))
var.clust <- hclust(as.dist(1-cor.mat))

#pdf("10x_CellLineMix_TestGroups.pdf")
#for (ngroups in 16:20) {
ngroups <- 20

# Correlation heatmap
hm1 <- Heatmap(cor.mat,
               cluster_columns = var.clust,
               cluster_rows = var.clust,
               row_split = ngroups,
               column_split = ngroups,
               row_gap = unit(0.5, "mm"),
               column_gap = unit(0.5, "mm"),
               row_names_gp = gpar(fontsize = 4),
               column_names_gp = gpar(fontsize = 4),
               row_title_gp = gpar(fontsize = 4),
               width = unit(100, "mm"),
               height = unit(100, "mm"),
               column_title = ngroups)
pdf(paste0("TenX_K562_clones_3_Correlation.pdf"))
print(hm1)
dev.off()

# Determine variants in each group as displayed by Heatmap package:
Groups.tib <- tibble(var = names(cutree(var.clust, k = ngroups)), Cut = cutree(var.clust, k = ngroups))[var.clust$order,]
Groups.tib <- Groups.tib %>% mutate(Group = match(Cut, unique(Cut)))
Groups.ls <- split(Groups.tib$var, Groups.tib$Group)
# Order from Groups with most cells to groups least cells
GroupSizes.ls <- lapply(Groups.ls, function(x) sum(af.dm[x,] > 1))
Groups.ls <- Groups.ls[rev(order(unlist(GroupSizes.ls)))]
# And then order from Groups with many to few variants
Groups.ls <- Groups.ls[order(-lengths(Groups.ls))]

### Customize column order
plot_order.mat <- af_subset.mat#[,order(-af_subset.mat["8213_G>A",])]
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
        tmp.cor[is.na(tmp.cor)] <- 1
        tmp.hclust <- hclust(as.dist(1-tmp.cor))
        tmp_order.mat <- plot_order.mat[,tmp.hclust$labels[tmp.hclust$order]]
        plot_order.mat <- cbind(tmp_order.mat, plot_order.mat[,!colnames(plot_order.mat) %in% colnames(tmp_order.mat)])
    }
}

# Generate a matrix to plot, while maintaining sorting order, adding back large variants, and restoring VAFs of 0-2%.
plot.mat <- af_subset.mat[rownames(plot_order.mat), colnames(plot_order.mat)]

# Plot
hm2 <- Heatmap(plot.mat,
               col = c("#FCFCFC","#FFEDB0","#FFDF5F","#FEC510","#FA8E24","#F14C2B","#DA2828","#BE2222","#A31D1D"),
               show_row_names = T,
               show_column_names = F,
               cluster_columns = F,
               cluster_rows = F,
               row_names_gp = gpar(fontsize = 5),
               name = "AF",
               heatmap_legend_param = list(at = c(0, 5, 10, 15), border = "#000000", grid_height = unit(10, "mm")),
               border = T,
               width = unit(100, "mm"),
               height = unit(100, "mm"),
               use_raster = T)
pdf(paste0("TenX_K562_clones_4_Heatmap.pdf"))
print(hm2)
#}
dev.off()

