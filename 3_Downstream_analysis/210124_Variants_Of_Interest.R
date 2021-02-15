# Peter van Galen, 210124
# Find a biologically interesting message in BPDCN712 sample


#~~~~~~~~~~~~~~~~~~~~~#
#### Prerequisites ####
#~~~~~~~~~~~~~~~~~~~~~#

options(stringsAsFactors = FALSE)
options(scipen = 999)

library(tidyverse)
library(data.table)
library(SummarizedExperiment)
library(Seurat)
library(Matrix)
library(ggrastr)
library(ComplexHeatmap)
library(gdata)
library(readxl)
library(ggrepel)
library(circlize)
library(gridExtra)
#library(GGally) # for ggpairs
#library(ggforce) # for geom_sina
#library(ggbeeswarm)
#library(stringr)

rm(list=ls())
setwd("~/DropboxPartners/Projects/Maester/AnalysisPeter/210123_BPDCN712_Diagnosis")

# Functions & colors
source("../201007_FunctionsGeneral.R")
popcol.df <- read.xls("~/DropboxPartners/Pipelines/AuxiliaryFiles/PopCol.xlsx", sheet = 3, row.names = 1)
mycol.ch <- popcol.df$hex
names(mycol.ch) <- rownames(popcol.df)
heatcol.ch <- read.xls("~/DropboxPartners/Pipelines/AuxiliaryFiles/PopCol.xlsx", sheet = 4, header = F)$V1

### Import data
# Load Seurat object
seu <- readRDS("BPDCN712_Seurat.rds")

# Load Maegtk, calculate allele frequencies
maegatk.rse <- readRDS("BPDCN712_Maegatk.rds")
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))*100

all(colnames(af.dm) == colnames(seu))

# Extract Seurat metadata
cells.tib <- as_tibble(seu@meta.data, rownames = "cell")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Which cell types to compare ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

cells.tib$CellType %>% table

CellSubsets.ls <- list(unionCells = cells.tib$cell,
                       HSPC = filter(cells.tib, CellType %in% c("HSC", "Prog"))$cell,
                       Ery = filter(cells.tib, CellType %in% c("EarlyE", "LateE"))$cell,
                       Myeloid = filter(cells.tib, CellType %in% c("ProMono", "Mono", "ncMono", "cDC", "pDC"))$cell,
                       Lymphoid = filter(cells.tib, CellType %in% c("T", "CTL"))$cell,
                       B = filter(cells.tib, CellType %in% c("B", "Plasma"))$cell,
                       NK = filter(cells.tib, CellType == "NK")$cell,
                       CTC = filter(cells.tib, BPDCN_Tumor > 0.75)$cell)
lengths(CellSubsets.ls)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Collect variant information ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# This can take a while so consider loading vars.tib from the bottom

# Get the mean allele frequency and coverage for every cell subset.
mean_af.ls <- lapply(CellSubsets.ls, function(x) rowMeans(af.dm[,x]))
mean_cov.ls <- lapply(CellSubsets.ls, function(x) rowMeans(assays(maegatk.rse)[["coverage"]][,x])[as.numeric(cutf(rownames(af.dm), d = "_"))])
names(mean_af.ls) <- paste0("mean_af.", names(mean_af.ls))
names(mean_cov.ls) <- paste0("mean_cov.", names(mean_cov.ls))

# Get the quantiles of the VAFs of each variant in each cell subset
quantiles <- c("q01" = 0.01, "q05" = 0.05, "q10" = 0.1, "q50" = 0.5,
               "q90" = 0.9, "q95" = 0.95, "q97" = 0.97, "q98" = 0.98, "q99" = 0.99, "q999" = 0.999)
# This can take a while
start_time <- Sys.time()
quantiles.ls <- lapply(quantiles, function(x) lapply(CellSubsets.ls, function(y) apply(af.dm[,y], 1, quantile, x) ))
Sys.time() - start_time

# Get the mean quality for each variant. This can take a few hours.
start_time <- Sys.time()
qual.num <- sapply(rownames(af.dm), function(x) {
    pos <- as.numeric( cutf(x, d = "_") )
    message(pos)
    mut <- cutf(x, d = ">", f = 2)
    # Only use cells in which the base was sequenced. Use reverse only because that's how we amplify transcripts.
    covered <- assays(maegatk.rse)[[str_c(mut, "_counts_rev")]][pos,] > 0
    # Get the mean quality for this call
    qual <- mean( assays(maegatk.rse)[[str_c(mut, "_qual_rev")]][pos, covered] )
    return(qual)
})
Sys.time() - start_time

# Collect all information in a tibble
vars.tib <- as_tibble(do.call(cbind, c(mean_af.ls, mean_cov.ls, unlist(quantiles.ls, recursive = F))), rownames = "var")
vars.tib <- add_column(vars.tib, quality = qual.num, .before = 2)

# Save for fast loading next time
write_tsv(vars.tib, "210124_all_variants.txt")
vars.tib <- read_tsv("210124_all_variants.txt")


#~~~~~~~~~~~~~~~~~~#
#### Background ####
#~~~~~~~~~~~~~~~~~~#

# Get mean allele frequence in K562 and BT142 cells
backgroundvars.tib <- read_tsv("../200922_Cell_line_mixes/TenX_CellLineMix_Variants1.txt")
backgroundvars.tib <- backgroundvars.tib %>% select(var, mean_cov.unionCells, mean_af.K562, mean_af.BT142)
vars.tib <- vars.tib %>% left_join(select(backgroundvars.tib, var, mean_af.K562, mean_af.BT142), by = "var")

# Function to determine the VAF fold difference between the two cell lines
fold_difference.var <- function(a, b) {
    # Distance to homoplasmy
    x <- min(a, abs(100-a))
    y <- min(b, abs(100-b))
    # Absolute difference
    d <- abs(a-b)
    # Fold difference
    z <- (min(x,y)+d) / min(x,y)
    return(z)
}

# Test if the function gives the right answer
fold_difference.var(80, 60)
fold_difference.var(2, 20)
fold_difference.var(1, NA)
fold_difference.var(1, 0.2)
fold_difference.var(99, 1)
fold_difference.var(99.9, 10)
fold_difference.var(99.9, 90)
fold_difference.var(50, 55)
fold_difference.var(45, 50)

backgroundvars.tib$FoldChange <- apply(backgroundvars.tib, 1, function(x) fold_difference.var(as.numeric(x[3]), as.numeric(x[4])))

blocklist.var <- backgroundvars.tib %>% filter(mean_cov.unionCells > 5) %>%
    filter(between(mean_af.K562, 0.1, 99.9), between(mean_af.BT142, 0.1, 99.9), FoldChange < 5) %>% .$var
blocklist.var <- c(blocklist.var, backgroundvars.tib[grepl("N", backgroundvars.tib$var),]$var)

# Polygon coordinates
tib <- tibble(x = c(0.1, 0.1, 80/4.8,  99.5, 99.9, 99.9, 5*100/6,  0.5, 0.1),
              y = c(0.1, 0.5, 5*100/6, 99.9, 99.9, 99.5, 80/4.8,   0.1, 0.1))

### To calculate where x and y intersect
# Upper left and lower right
# line 1: x = 0, y = 0, slope = 5
# line 2: x = 0, y = 80, slope = 0.2
# y = mx + b (m = slope, b = y intersect)
# y = 5x + 0 (slope of the first line)
# y = 0.2x + 80
# 5x + 0 = 0.2x + 80
# 5x - 0.2x = 80
# 4.8x = 80
# x = 80/4.8

# Close to the corners
# slope = (y2-y1) / (x2-x1)
# 0.2 = (99.9-80) / (x2-0)
# 0.2 = 19.9 / x2
# 0.2*x = 19.9
# x = 19.9/0.2
# 99.5

points.tib <- backgroundvars.tib %>% filter(mean_cov.unionCells > 5) %>%
    mutate(block = var %in% blocklist.var) %>% arrange(block)

pdf("210124_1_Blocklist.pdf")
for (lim in list(c(0,100), c(0,1), c(99,100))) {
    print(
        ggplot(data = tib, mapping = aes(x = x, y = y)) +
            geom_polygon(fill = "#dc143c33") +
            geom_hline(yintercept = c(0.1,99.9), col = "#778899") +
            geom_vline(xintercept = c(0.1, 99.9), col = "#778899") +
            geom_abline(slope = c(0.2,5), intercept = 0, col = "#778899") +
            geom_abline(slope = c(0.2,5), intercept = c(80,-400), col = "#778899") +
            geom_point(data = points.tib, mapping = aes(x = mean_af.K562, y = mean_af.BT142, color = block)) +
            scale_color_manual(values = c("#111111", "#dc143c")) +
            coord_cartesian(xlim = unlist(lim), ylim = unlist(lim)) +
            ylab("Mean VAF in K562") +
            xlab("Mean VAF in BT142") +
            theme(aspect.ratio = 1)
    )
}
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Look into variants ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~#

### Small expansions with high VAF at diagnosis (see also 201015_Variants_Of_Interest.R)
save_name <- "210124_Dx"
voi.ch <- vars.tib %>% filter(mean_cov.unionCells > 5, quality >= 30, ! var %in% blocklist.var) %>% 
    filter(q90.unionCells == 0, q999.unionCells > 50) %>% .$var
# This one is background noise
voi.ch <- voi.ch[!grepl("1583_A>G", voi.ch)]
# Transitions
str_view(voi.ch, "G>A|A>G|C>T|T>C"); mean( str_count(voi.ch, "G>A|A>G|C>T|T>C") )
# Compare VAF between lymphoid and myeloid
vars.tib %>% filter(mean_cov.unionCells > 5, quality >= 30, ! var %in% blocklist.var) %>% arrange(var %in% voi.ch) %>%
    mutate(difference = abs(mean_af.Myeloid-mean_af.Lymphoid)) %>%
    ggplot(aes(x = mean_af.Myeloid, y = mean_af.Lymphoid, color = var %in% voi.ch,
               label = ifelse(var %in% voi.ch & difference > 2, var, ""))) +
    geom_point() +
    geom_text_repel(show.legend = F) +
    coord_cartesian(xlim = c(0,15), ylim = c(0,15)) +
    labs(color = "Condition") +
    theme(aspect.ratio = 1)

### Variants that might distinguish CTCs from non-malignant bone marrow cells
save_name <- "210124_CTCs"
voi.ch <- vars.tib %>% filter(mean_cov.unionCells > 5, quality >= 30, ! var %in% blocklist.var) %>%
    filter(mean_af.CTC > 0.5, mean_af.unionCells < 0.1) %>% .$var
# Transitions
str_view(voi.ch, "G>A|A>G|C>T|T>C"); mean( str_count(voi.ch, "G>A|A>G|C>T|T>C") )
vars.tib %>% filter(mean_cov.unionCells > 5, quality >= 30, ! var %in% blocklist.var) %>% arrange(var %in% voi.ch) %>%
    ggplot(aes(x = mean_af.unionCells, y = mean_af.CTC, color = var %in% voi.ch, label = ifelse(var %in% voi.ch, var, ""))) +
    geom_point() +
    geom_label_repel(show.legend = F) +
    coord_cartesian(xlim = c(0,5), ylim = c(0,5)) +
    labs(color = "Condition") +
    theme(aspect.ratio = 1)
# As you can see in the plots below, these "CTC" variants are mostly likely pure garbage. Each of them only occurs in one pDC and also in bunch of other cells, i.e., they could just be rare false positives.


#~~~~~~~~~~~~~~~~~~#
#### Plot UMAPs ####
#~~~~~~~~~~~~~~~~~~#

pdf(str_c(save_name, "_2_UMAPs.pdf"))
for (v in voi.ch) {
    message(v)
    
    # Add info of variant of interest
    cells.tib$af_voi <- af.dm[v,cells.tib$cell]
    cells.tib$cov_voi <- assays(maegatk.rse)[["coverage"]][as.numeric( cutf(v, d = "_") ),cells.tib$cell]
    
    print(
        cells.tib %>% arrange(af_voi) %>% filter(cov_voi > 5) %>%
            ggplot(aes(x = UMAP_1, y = UMAP_2, color = af_voi)) + # change to cov_voi to see coverage
            geom_point_rast(size = 1) +
            scale_color_gradientn(colors = heatcol.ch[4:10], limits = c(0,100), n.breaks = 3) +
            theme_classic() +
            theme(aspect.ratio = 1, axis.line = element_blank(), plot.title = element_text(hjust=0.5),
                  panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
            ggtitle(v)
    )
    print(
        cells.tib %>% arrange(af_voi) %>% mutate(Rank_sorted_cells = row_number()) %>%
            ggplot(aes(x = Rank_sorted_cells, y = af_voi, color = CellType)) +
            scale_color_manual(values = popcol.df[levels(cells.tib$CellType), "hex"]) +
            geom_point_rast() + 
            ylab("Variant allele frequency") +
            theme(aspect.ratio = 0.5, plot.title = element_text(hjust=0.5)) +
            ggtitle(v)
    )
    cells.tib$af_voi <- NULL
    cells.tib$cov_voi <- NULL
}
dev.off()

# See below UMAP plots for supplementary data


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Group similar variants ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

### Prepare matrix of variants of interest in cells that are positive for at least one
af_voi.mat <- af.dm[voi.ch,]
af_subset.mat <- af_voi.mat[,apply(af_voi.mat, 2, function(x) sum(x > 1) > 0)]
dim(af.dm); dim(af_voi.mat); dim(af_subset.mat)
af_subset.mat[,1:2]

### Variant correlation & cluster
cor.mat <- cor(t(af_subset.mat))
var.clust <- hclust(as.dist(1-cor.mat))

# Correlation heatmap
plot(var.clust$height, ylim = c(0, max(var.clust$height))); abline(h = 0.5, col = "red")
ngroups <- sum(var.clust$height > 0.5) + 1 # not sure how consistently this performs but it works for small_dx_clones

hm1 <- Heatmap(cor.mat,
               col = colorRamp2(c(-1,0,1), c("blue", "#DDDDDD", "red")),
               cluster_columns = var.clust,
               cluster_rows = var.clust,
               row_split = switch(ngroups < length(voi.ch), ngroups),
               column_split = switch(ngroups < length(voi.ch), ngroups),
               show_row_dend = F, # without this the visualizationn does not complete
               show_column_dend = F, # without this the visualizationn does not complete
               row_gap = unit(0.5, "mm"),
               column_gap = unit(0.5, "mm"),
               row_names_gp = gpar(fontsize = 10),
               column_names_gp = gpar(fontsize = 10),
               row_title_gp = gpar(fontsize = 10),
               width = unit(100, "mm"),
               height = unit(100, "mm"),
               column_title = ngroups)
pdf(str_c(save_name, "_3_Correlation.pdf"))
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

# Reorder voi.ch (after making sure none will be lost):
identical(sort(voi.ch), sort(unlist(str_split(Groups.tib$vars, ", "))))
voi.ch <- unlist(str_split(Groups.tib$vars, ", "))


#~~~~~~~~~~~~~~~~~~~#
#### VAF heatmap ####
#~~~~~~~~~~~~~~~~~~~#

### Customize column order
# Sort for all variants from the correlation matrix
plot_order.mat <- af_subset.mat[unlist(str_split(Groups.tib$vars, ", ")),]
# OPTION 1 (better for K562 clones)
# Only use VAFs of >2% for sorting
#plot_order.mat[plot_order.mat < 2 ] <- 0
# Order from high to low
#for (x in rev(strsplit(Groups.tib$vars, ", "))) {
#    if (length(x) == 1) { # for individual variants, order by VAF
#        plot_order.mat <- plot_order.mat[,order(-plot_order.mat[x,])]
#    } else { # for variants in a Group, use Pearson clustering for ordering
#        tmp.mat <- plot_order.mat[x,! colSums(plot_order.mat[x,]) == 0]
#        tmp.cor <- cor(tmp.mat, method = "pearson")
#        tmp.cor[is.na(tmp.cor)] <- 1
#        tmp.hclust <- hclust(as.dist(1-tmp.cor))
#        tmp_order.mat <- plot_order.mat[,tmp.hclust$labels[tmp.hclust$order]]
#        plot_order.mat <- cbind(tmp_order.mat, plot_order.mat[,!colnames(plot_order.mat) %in% colnames(tmp_order.mat)])
#    }
#}
# OPTION 2 (better for Dx clones)
for (x in rev(strsplit(Groups.tib$vars, ", "))) {
    if (length(x) == 1) {
        plot_order.mat <- plot_order.mat[,order(-plot_order.mat[x,])]
    } else {
        plot_order.mat <- plot_order.mat[,order(-colSums(plot_order.mat[x,]))]
    }
}

# Generate a matrix to plot, while maintaining sorting order, adding back large variants, and restoring VAFs of 0-2%.
plot.mat <- af_subset.mat[rownames(plot_order.mat), colnames(plot_order.mat)]

# Add annotation bars
anno.tib <- tibble(cell = colnames(plot.mat)) %>% left_join(cells.tib, by = "cell") %>%
    select(CellType, ASXL1.G642fs, `TET2.S792*`, `TET2.Q1034*`, `TET2.R1216*`, TET2.H1380Y)
colors.ch <- popcol.df$hex
names(colors.ch) <- rownames(popcol.df)
ha <- HeatmapAnnotation(df = data.frame(anno.tib),
                        col = list(CellType = colors.ch,
                                   ASXL1.G642fs = c(wildtype = "#FFFFFF", mutant = "#000000", `no call` = "#a0a0a0"),
                                   TET2.S792. = c(wildtype = "#FFFFFF", mutant = "#000000", `no call` = "#a0a0a0"),
                                   TET2.Q1034. = c(wildtype = "#FFFFFF", mutant = "#000000", `no call` = "#a0a0a0"),
                                   TET2.R1216. = c(wildtype = "#FFFFFF", mutant = "#000000", `no call` = "#a0a0a0"),
                                   TET2.H1380Y = c(wildtype = "#FFFFFF", mutant = "#000000", `no call` = "#a0a0a0")))
# Take out top 1% of values
#plot.mat[plot.mat>quantile(plot.mat, 0.99)] <- quantile(plot.mat, 0.99)

# Plot
hm2 <- Heatmap(plot.mat,
               col = colorRamp2(seq(0, round(max(plot.mat)), length.out = 9),
                                c("#FCFCFC","#FFEDB0","#FFDF5F","#FEC510","#FA8E24","#F14C2B","#DA2828","#BE2222","#A31D1D")),
               show_row_names = T,
               show_column_names = F,
               cluster_columns = F,
               cluster_rows = F,
               row_names_gp = gpar(fontsize = 10),
               name = "AF",
               heatmap_legend_param = list(#at = seq(0, round(max(plot.mat)), length.out = 5),
                   border = "#000000", grid_height = unit(10, "mm")),
               top_annotation = ha,
               border = T,
               #width = unit(100, "mm"),
               #height = unit(100, "mm"),
               use_raster = T)
pdf(str_c(save_name, "_4_Heatmap.pdf"), width = 12, height = 6)
print(hm2)
dev.off()

# Save variants
write_tsv(tibble(var = voi.ch), file = str_c(save_name, "_vois.txt"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### UMAPs for supplement ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Save UMAPs: all on one page for supplemental figure.
# Make everything 1.5x the final size and then reduce size in Illustrator. This was necessary because the size argument in geom_point didn't work properly.
p <- list()
for (v in setdiff(voi.ch, c("683_G>A", "2593_G>A"))) {
    #v <- "4104_A>G"
    cells.tib$af_voi <- af.dm[v,cells.tib$cell]
    cells.tib$cov_voi <- assays(maegatk.rse)[["coverage"]][as.numeric( cutf(v, d = "_") ),cells.tib$cell]
    
    p[[v]] <- cells.tib %>% arrange(af_voi) %>% #filter(cov_voi > 5) %>%
        ggplot(aes(x = UMAP_1, y = UMAP_2, color = af_voi)) +
        geom_point_rast(size = 0.1, raster.dpi = 900) +
        scale_color_gradientn(colors = heatcol.ch[4:10], limits = c(0,100), n.breaks = 3) +
        theme_classic() +
        theme(aspect.ratio = 1, axis.line = element_blank(), plot.title = element_text(hjust=0.5, size = 9),
              panel.border = element_rect(colour = "black", fill=NA, size=0.5/2.13*3/2),
              legend.position = "none", axis.title = element_blank(), axis.text = element_blank(),
              axis.ticks = element_blank()) +
        ggtitle(v)
    
    cells.tib$af_voi <- NULL
    cells.tib$cov_voi <- NULL
}

pdf(file = str_c(save_name, "_2_UMAPs_combined.pdf"), height = (11-2)*1.5, width = (8.5-2)*1.5)
grid.arrange(grobs = p, ncol = 4)
dev.off()
