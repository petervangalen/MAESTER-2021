# Peter van Galen, 210124
# Assess clonal relationships in the clonal hematopoiesis sample


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
library(readxl)
library(ggrepel)
library(circlize) # for colorRamp2

rm(list=ls())
setwd("~/DropboxPartners/Projects/Maester/AnalysisPeter/4_CH_sample")

# Functions and colors (available at https://github.com/vangalenlab/MAESTER-2021)
source("../210215_FunctionsGeneral.R")
popcol.df <- read_excel("../MAESTER_colors.xlsx")
mycol.ch <- popcol.df$hex
names(mycol.ch) <- popcol.df$name
heatcol.ch <- read_excel("../MAESTER_colors.xlsx", sheet = 2, col_names = "heatcol")$heatcol

# Load Seurat object (available at https://vangalenlab.bwh.harvard.edu/maester-2021/)
seu <- readRDS("BPDCN712_Seurat.rds")

# Load Maegtk, calculate allele frequencies (https://vangalenlab.bwh.harvard.edu/maester-2021/))
maegatk.rse <- readRDS("BPDCN712_Maegatk.rds")
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))*100
# Check (should be TRUE)
all(colnames(af.dm) == colnames(seu))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Which cell types to compare ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Extract Seurat metadata
cells.tib <- as_tibble(seu@meta.data, rownames = "cell")
cells.tib$CellType %>% table

CellSubsets.ls <- list(unionCells = cells.tib$cell,
                       HSPC = filter(cells.tib, CellType %in% c("HSC", "Prog"))$cell,
                       Ery = filter(cells.tib, CellType %in% c("EarlyE", "LateE"))$cell,
                       Myeloid = filter(cells.tib, CellType %in% c("ProMono", "Mono", "ncMono", "cDC", "pDC"))$cell,
                       Lymphoid = filter(cells.tib, CellType %in% c("T", "CTL"))$cell,
                       B = filter(cells.tib, CellType %in% c("B", "Plasma"))$cell,
                       NK = filter(cells.tib, CellType == "NK")$cell,
                       CTC = filter(cells.tib, BPDCN_Tumor_score > 0.75)$cell)
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
backgroundvars.tib <- read_tsv("../3_Cell_line_mixes_variants/TenX_CellLineMix_Variants1.txt")
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

# Parameters to select small clones with high VAF
voi.ch <- vars.tib %>% filter(mean_cov.unionCells > 5, quality >= 30, ! var %in% blocklist.var) %>% 
    filter(q90.unionCells == 0, q999.unionCells > 50) %>% .$var
# This one is background noise
voi.ch <- voi.ch[!grepl("1583_A>G", voi.ch)]
# Assess transitions vs. transversions
str_view(voi.ch, "G>A|A>G|C>T|T>C"); mean( str_count(voi.ch, "G>A|A>G|C>T|T>C") )

# Compare VAF between lymphoid and myeloid
vars.tib %>% filter(mean_cov.unionCells > 5, quality >= 30, ! var %in% blocklist.var) %>% arrange(var %in% voi.ch) %>%
    mutate(difference = abs(mean_af.Myeloid-mean_af.Lymphoid)) %>%
    ggplot(aes(x = mean_af.Myeloid, y = mean_af.Lymphoid, color = var %in% voi.ch,
               label = ifelse(var %in% voi.ch & difference > 2, var, ""))) +
    geom_point() +
    geom_text_repel(show.legend = F) +
    coord_cartesian(xlim = c(0,15), ylim = c(0,15)) +
    labs(color = "Variant of interest") +
    theme_classic() +
    theme(aspect.ratio = 1)


#~~~~~~~~~~~~~~~~~~#
#### Plot UMAPs ####
#~~~~~~~~~~~~~~~~~~#

pdf(str_c("210124_2_UMAPs.pdf"))
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
            scale_color_manual(values = mycol.ch) +
            geom_point_rast() + 
            ylab("Variant allele frequency") +
            theme(aspect.ratio = 0.5, plot.title = element_text(hjust=0.5)) +
            ggtitle(v)
    )
    cells.tib$af_voi <- NULL
    cells.tib$cov_voi <- NULL
}
dev.off()

# See below for UMAP plots for supplementary data


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

# Assess how correlated variants are and group them together
plot(var.clust$height, ylim = c(0, max(var.clust$height)))
# Make 23 groups of variants, i.e. 20 alone and 3 groups of two variants. This is determined empirically.
ngroups <- 23

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
pdf(str_c("210124_3_Correlation.pdf"))
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

# Sort for all variants from the correlation matrix
plot_order.mat <- af_subset.mat[unlist(str_split(Groups.tib$vars, ", ")),]

# Customize column order. This is different from the strategy for K562 subclones.
for (x in rev(strsplit(Groups.tib$vars, ", "))) {
    if (length(x) == 1) {
        plot_order.mat <- plot_order.mat[,order(-plot_order.mat[x,])]
    } else {
        plot_order.mat <- plot_order.mat[,order(-colSums(plot_order.mat[x,]))]
    }
}

# Add annotation bars
anno.tib <- tibble(cell = colnames(plot_order.mat)) %>% left_join(cells.tib, by = "cell") %>%
    select(CellType, ASXL1.G642fs, `TET2.S792*`, `TET2.Q1034*`, `TET2.R1216*`, TET2.H1380Y)
ha <- HeatmapAnnotation(df = data.frame(anno.tib),
                        col = list(CellType = mycol.ch,
                                   ASXL1.G642fs = c(wildtype = "#FFFFFF", mutant = "#000000", `no call` = "#a0a0a0"),
                                   TET2.S792. = c(wildtype = "#FFFFFF", mutant = "#000000", `no call` = "#a0a0a0"),
                                   TET2.Q1034. = c(wildtype = "#FFFFFF", mutant = "#000000", `no call` = "#a0a0a0"),
                                   TET2.R1216. = c(wildtype = "#FFFFFF", mutant = "#000000", `no call` = "#a0a0a0"),
                                   TET2.H1380Y = c(wildtype = "#FFFFFF", mutant = "#000000", `no call` = "#a0a0a0")))

# Plot
hm2 <- Heatmap(plot_order.mat,
               col = colorRamp2(seq(0, round(max(plot_order.mat)), length.out = 9),
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
pdf(str_c("210124_4_Heatmap.pdf"), width = 12, height = 6)
print(hm2)
dev.off()

# Save variants
write_tsv(tibble(var = voi.ch), file = str_c("210124_vois.txt"))


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

pdf(file = str_c("210124_2_UMAPs_combined.pdf"), height = (11-2)*1.5, width = (8.5-2)*1.5)
grid.arrange(grobs = p, ncol = 4)
dev.off()


#~~~~~~~~~~~~~~~~~~#
#### More plots ####
#~~~~~~~~~~~~~~~~~~#

# Mean depth in each cell type
mean_depth.tib <- tibble(cell = colnames(maegatk.rse),
                         depth = maegatk.rse@colData$depth)
pdf("210124_5_Coverage_per_CellType.pdf")
mean_depth.tib %>% left_join(select(cells.tib, cell, CellType), by = "cell") %>%
    group_by(CellType) %>% summarize(mean_depth = mean(depth)) %>%
    ggplot(aes(x = CellType, y = mean_depth)) +
    geom_bar(stat = "identity") +
    ylab("Mean depth") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
dev.off()


# This can be assessed for the different variants as follows
subset_cov.dgm <- assays(maegatk.rse)[["coverage"]][as.numeric(cutf(voi.ch, d = "_")),]
rownames(subset_cov.dgm) <- voi.ch
cov.tib <- as_tibble(t(as.matrix(subset_cov.dgm)), rownames = "cell") %>%
    left_join(select(cells.tib, cell, CellType), by = "cell") %>% relocate(CellType, .before = 2)
pdf("210124_6_Coverage_per_Variant.pdf", width = 12)
cov.tib %>% pivot_longer(cols = -c(cell, CellType), names_to = "variant", values_to = "cov") %>%
    group_by(variant, CellType) %>% summarize(mean_cov = mean(cov)) %>%
    mutate(variant = factor(variant, levels = voi.ch)) %>%
    ggplot(aes(x = CellType, y = mean_cov)) +
    geom_bar(stat = "identity") +
    facet_wrap(~variant, scales = "free_y") +
    ylab("Mean coverage") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
dev.off()


# Similarly, this information can be extracted from vars.tib for each cell subset
vars.tib %>% select(var, contains("mean_cov")) %>% .[match(voi.ch, vars.tib$var),] %>%
    mutate(var = factor(var, levels = var)) %>%
    pivot_longer(cols = -var, names_to = "CellType", values_to = "mean_coverage") %>%
    mutate(CellType = gsub("mean_cov.", "", CellType)) %>%
    filter(! CellType %in% c("unionCells", "CTC")) %>%
    mutate(CellType = factor(CellType, levels = unique(CellType))) %>%
    ggplot(aes(x = CellType, y = mean_coverage)) +
    geom_bar(stat = "identity") +
    facet_wrap(~var, scales = "free_y") +
    ylab("Mean coverage") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())


### Plot of variant locations
voi_order.ch <- voi.ch[order(as.numeric(cutf(voi.ch, d = "_")))]
voi_along_chrM.num <- rep(0, 16569)
voi_along_chrM.num[as.numeric(cutf(voi_order.ch, d = "_"))] <- 1
labels.ch <- rep("", 16569)
labels.ch[as.numeric(cutf(voi_order.ch, d = "_"))] <- voi_order.ch

# Colors
varcol.ch <- mycol.ch[sapply(voi_order.ch, function(x) grep(x, names(mycol.ch)))]
data.frame(names(varcol.ch), voi_order.ch)
names(varcol.ch) <- voi_order.ch
barcol.ch <- varcol.ch[match(labels.ch, names(varcol.ch))]

# Genes
GenePos.tib <- tibble(Names = c("MT.ATP6", "MT.ATP8", "MT.COX1", "MT.COX2", "MT.COX3", "MT.CYTB", "MT.ND1", "MT.ND2", "MT.ND3",
                                "MT.ND4", "MT.ND4L", "MT.ND5", "MT.ND6", "MT.RNR1", "MT.RNR2"),
                      start = c(8527, 8366, 5904, 7586, 9207, 14747, 3307, 4470, 10059, 10760, 10470, 12337, 14149, 648, 1671), 
                      end = c(9207, 8572, 7445, 8269, 9990, 15887, 4262, 5511, 10404, 12137, 10766, 14148, 14673, 1601, 3229))
GenePos.tib <- GenePos.tib %>% arrange(start) %>%
    mutate(mid = round((end-start)/2+start,0), ycoord = rep(c(-0.05,-0.07), length.out = 15))

pdf("210124_7_Variant_locations.pdf", height = 3, width = 10)
par(xpd=T)
barplot(voi_along_chrM.num, border = barcol.ch,
        yaxt = "n", xlab = "position", cex.names = 0, space = 0, ylim = c(0,1.5))
axis(side = 1, at = c(0, 16569))
points(x = as.numeric(cutf(voi_order.ch, d = "_")), y = rep(1,length(voi_order.ch)),
       pch = 16, col = varcol.ch)

for (n in 1:nrow(GenePos.tib)) {
    lines(x = c(GenePos.tib$start[n], GenePos.tib$end[n]), y = c(GenePos.tib$ycoord[n], GenePos.tib$ycoord[n]), col = "blue")
    text(x = GenePos.tib$mid[n], y = GenePos.tib$ycoord[n],
         labels = gsub("MT.", "", GenePos.tib$Names[n]), pos = 1, col = "blue", cex = 0.7)
}

textcol.ch <- na.omit(barcol.ch)
attributes(textcol.ch) <- NULL
text(x = as.numeric(cutf(voi_order.ch, d = "_")), y = 1.1, labels = voi_order.ch, pos = 4, srt = 45,
     col = varcol.ch, cex = 0.7)
dev.off()



