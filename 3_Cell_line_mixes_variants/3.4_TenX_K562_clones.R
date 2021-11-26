# Peter van Galen, 200130
# Detect subclones in K562 cells from the 10X cell line mixing experiment
# This is similar to 3.2_SW_K562_clones.R

#~~~~~~~~~~~~~~~~~~~~~#
#### Prerequisites ####
#~~~~~~~~~~~~~~~~~~~~~#

options(scipen = 999)

library(tidyverse)
library(SummarizedExperiment)
library(Seurat)
library(data.table)
library(Matrix)
library(ComplexHeatmap)
library(readxl)
library(ggrastr)

setwd("~/DropboxMGB/Projects/Maester/AnalysisPeter/3_Cell_line_mixes_variants")
rm(list=ls())

# Functions and colors (available at https://github.com/vangalenlab/MAESTER-2021)
source("../210215_FunctionsGeneral.R")
popcol.df <- read_excel("../MAESTER_colors.xlsx")
heatcol.ch <- read_excel("../MAESTER_colors.xlsx", sheet = 2, col_names = "heatcol")$heatcol

# Import MAEGATK data for 10X Cell Line Mix (available at https://vangalenlab.bwh.harvard.edu/maester-2021/)
maegatk.rse <- readRDS("../1_MT_Coverage/TenX_CellLineMix_All_mr3_maegatk.rds")

# This tibble contains UMAP coordinates and cell type classifications
# Available at https://github.com/vangalenlab/MAESTER-2021; generated in 3.3_TenX_CellLineMix_variants.R
cells.tib <- read_rds("TenX_CellLineMix_Cells.rds")

    
#~~~~~~~~~~~~~~~~~~~~#
#### Prepare data ####
#~~~~~~~~~~~~~~~~~~~~#

# Put maegatk object in the same order as cells.tib
maegatk.rse <- maegatk.rse[,cells.tib$cell]

# Prepare allele frequency matrix
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))*100

# Group cell IDs by cell type classification. Use the cell type classification from maegatk, not RNA-seq.
table( as.character(cells.tib$CellType_RNA) == as.character(cells.tib$CellType_MT) )
CellSubsets.ls <- list(unionCells = cells.tib$cell,
                       K562 = filter(cells.tib, CellType_MT == "K562")$cell,
                       BT142 = filter(cells.tib, CellType_MT == "BT142")$cell)

# Check coverage between cell types
cells.tib %>% mutate(coverage = maegatk.rse$depth) %>%
    ggplot(aes(x = CellType_RNA, y = coverage)) +
    geom_violin() +
    theme(aspect.ratio = 1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Collect variant information ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# This section takes a while; consider skipping to read_tsv at the end.
start_time <- Sys.time()

# Get the mean allele frequency and coverage for every cell subset.
mean_af.ls <- lapply(CellSubsets.ls, function(x) rowMeans(af.dm[,x]))
mean_cov.ls <- lapply(CellSubsets.ls, function(x) rowMeans(assays(maegatk.rse)[["coverage"]][,x])[as.numeric(cutf(rownames(af.dm), d = "_"))])
names(mean_af.ls) <- paste0("mean_af.", names(mean_af.ls))
names(mean_cov.ls) <- paste0("mean_cov.", names(mean_cov.ls))

# Get the quantiles of the VAFs of each variant in each cell subset. This takes 2.5 minutes.
quantiles <- c("q01" = 0.01, "q10" = 0.1, "q50" = 0.5, "q90" = 0.9, "q99" = 0.99)
quantiles.ls <- lapply(quantiles, function(x) lapply(CellSubsets.ls, function(y) apply(af.dm[,y], 1, quantile, x) ))
Sys.time() - start_time

# Get the mean quality for each variant.
assays.ls <- lapply(maegatk.rse@assays$data, function(x) as.matrix(x))
start_time <- Sys.time()
qual.num <- sapply(rownames(af.dm), function(x) {
    #x <- "2141_T>C"
    pos <- as.numeric( cutf(x, d = "_") )
    mut <- cutf(x, d = ">", f = 2)
    # Get the mean quality of reads for this call (only use cells in which the base was sequenced) - forward
    covered_fw <- assays.ls[[str_c(mut, "_counts_fw")]][pos,] > 0
    qual_fw <- assays.ls[[str_c(mut, "_qual_fw")]][pos, covered_fw]
    # Same for reverse
    covered_rev <- assays.ls[[str_c(mut, "_counts_rev")]][pos,] > 0
    qual_rev <- assays.ls[[str_c(mut, "_qual_rev")]][pos, covered_rev]
    qual <- mean(c(qual_fw, qual_rev))
    return(qual)
})
Sys.time() - start_time

# Collect all information in a tibble
vars.tib <- as_tibble(do.call(cbind, c(mean_af.ls, mean_cov.ls, unlist(quantiles.ls, recursive = F))), rownames = "var")
vars.tib <- add_column(vars.tib, quality = qual.num, .before = 2)

# Plot quality metrics
vars.tib %>% arrange(quality) %>% mutate(key = row_number()) %>%
    mutate(mut = cutf(var, d = ">", f = 2)) %>%
    ggplot(aes(x = key, y = quality, color = mut)) +
    geom_point() +
    xlab(label = "Rank sorted variant") +
    theme(aspect.ratio = 0.5)

# Save for fast loading next time
write_tsv(vars.tib, "TenX_CellLineMix_Variants2.txt")
vars.tib <- read_tsv("TenX_CellLineMix_Variants2.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Detect variants of interest ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

### Get the most interesting subclonal variants
voi.ch <- vars.tib %>%
    filter(mean_cov.unionCells > 100, q01.K562 < 1, q99.K562 > 2,
           mean_af.K562 > (5*mean_af.BT142),
           quality > 30) %>%
    .$var
# Add 2141_T>C as a positive control / homoplasmic variant
voi.ch <- c(voi.ch, "2141_T>C")
# The following variants were not detected in orthogonal bulk ATAC-seq analysis, so I will remove them
voi.ch <- setdiff(voi.ch, c("7693_C>T", "8251_G>A"))
# Assess transitions vs. transversions
str_view(voi.ch, "G>A|A>G|C>T|T>C"); mean( str_count(voi.ch, "G>A|A>G|C>T|T>C") )
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

### Prepare matrix of variants of interest in K562 cells with decent coverage
cells.tib$n_covered <- cells.tib %>% select(contains("cov_")) %>% apply(1, function(x) sum(x > 3))
af_subset.mat <- af.dm[voi.ch,filter(cells.tib, CellType_MT == "K562", n_covered == length(voi.ch))$cell]
# Cap at 15% VAF
af_subset.mat[af_subset.mat > 15] <- 15

### Variant correlation.
var_cor.mat <- af_subset.mat[! rownames(af_subset.mat) == "2141_T>C",]
cor.mat <- cor(t(var_cor.mat))
# Cluster
var.clust <- hclust(as.dist(1-cor.mat))

# Assess how correlated variants are and group them together
plot(var.clust$height, ylim = c(0, max(var.clust$height)))
# Make groups of variants. This number is determined empirically.
ngroups <- 15

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
Groups.tib <- Groups.tib %>% group_by(Group) %>% summarize(vars = toString(var), nvar = n())
GroupIDs.ls <- lapply(str_split(Groups.tib$vars, ", "), function(x) c(sapply(x, function(y) colnames(af.dm[,af.dm[y,] > 1]))))
Groups.tib$ncells <- unlist(lapply(GroupIDs.ls, function(x) length(unique(unlist(x)))))
# Order from Groups with most cells to groups with least cells, then from Groups with many to few variants.
Groups.tib <- Groups.tib %>% arrange(desc(ncells), desc(nvar))

### Customize column order
plot_order.mat <- af_subset.mat[unlist(str_split(Groups.tib$vars, ", ")),]
# Only use VAFs of >2% for sorting
plot_order.mat[plot_order.mat < 2 ] <- 0
# Order columns
for (x in rev(strsplit(Groups.tib$vars, ", "))) {
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

# Generate a matrix to plot, while maintaining sorting order, adding back 2141_T>C, and restoring VAFs of 0-2%.
plot.mat <- af_subset.mat[c("2141_T>C", rownames(plot_order.mat)), colnames(plot_order.mat)]

# Remove uninformative cells
plot_reduced.mat <- plot.mat[,colSums(plot.mat[! rownames(plot.mat) == "2141_T>C",]) > 2]
dim(plot.mat); dim(plot_reduced.mat)

# Plot
hm2 <- Heatmap(plot_reduced.mat,
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
dev.off()

# Save data for comparison with bulk ATAC-seq
vars.tib %>% filter(var %in% voi.ch) %>%
    select(var, quality, mean_af.unionCells, mean_af.K562, mean_af.BT142, mean_cov.unionCells, mean_cov.K562, mean_cov.BT142) %>%
    write_tsv(file = "FigureS7E_vars.txt")


