# Peter van Galen, 211113
# Assess clonal relationships in the clonal hematopoiesis sample: generate Figure 2B, Supplemental Figure 12


#~~~~~~~~~~~~~~~~~~~~~#
#### Prerequisites ####
#~~~~~~~~~~~~~~~~~~~~~#

options(scipen = 999)

library(tidyverse)
library(data.table)
library(SummarizedExperiment)
library(Seurat)
library(Matrix)
library(ggrastr)
library(ComplexHeatmap)
library(circlize) # for colorRamp2
library(readxl)
library(ggrepel)
library(gridExtra)

rm(list=ls())
setwd("~/DropboxMGB/Projects/Maester/AnalysisPeter/4_CH_sample")

# Functions and colors (available at https://github.com/vangalenlab/MAESTER-2021)
source("../Auxiliary_files/210215_FunctionsGeneral.R")
popcol.df <- read_excel("../MAESTER_colors.xlsx")
mycol.ch <- popcol.df$hex
names(mycol.ch) <- popcol.df$name
heatcol.ch <- read_excel("../MAESTER_colors.xlsx", sheet = 2, col_names = "heatcol")$heatcol

# Load Seurat object (available at https://vangalenlab.bwh.harvard.edu/maester-2021/)
seu <- readRDS("BPDCN712_Seurat_Final.rds")
metadata.tib <- as_tibble(seu@meta.data, rownames = "cell")

# Load Maegtk, calculate allele frequencies (https://vangalenlab.bwh.harvard.edu/maester-2021/))
maegatk.rse <- readRDS("BPDCN712_Maegatk_Final.rds")
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))*100
# Check (should be TRUE)
all(colnames(af.dm) == colnames(seu))

# Load 26 variants of interest from 4.2_Variant_Selection.R
voi.ch <- read.table("4.2_vois.txt")$V1

#~~~~~~~~~~~~~~~~~~~~~~~~#
# Group similar variants #
#~~~~~~~~~~~~~~~~~~~~~~~~#

# Some variants show similar clonal structures. Here I group these similar variants into clones.
# Prepare matrix of variants of interest in cells that are positive for at least one
af_voi.mat <- af.dm[voi.ch,]
af_subset.mat <- af_voi.mat[,apply(af_voi.mat, 2, function(x) sum(x > 1) > 0)]
dim(af.dm); dim(af_voi.mat); dim(af_subset.mat)

# Variant correlation & cluster
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
pdf(str_c("4.3_1_Correlation.pdf"))
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


#~~~~~~~~~~~~~#
# VAF heatmap #
#~~~~~~~~~~~~~#

# Sort for all variants from the correlation matrix
plot_order.mat <- af_subset.mat[unlist(str_split(Groups.tib$vars, ", ")),]

# Customize column order.
for (x in rev(strsplit(Groups.tib$vars, ", "))) {
    if (length(x) == 1) {
        plot_order.mat <- plot_order.mat[,order(-plot_order.mat[x,])]
    } else {
        plot_order.mat <- plot_order.mat[,order(-colSums(plot_order.mat[x,]))]
    }
}

# Add annotation bars (ASXL1 and TET2 mutations are removed from the final figure)
anno.tib <- tibble(cell = colnames(plot_order.mat)) %>% left_join(metadata.tib, by = "cell") %>%
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
               heatmap_legend_param = list(border = "#000000", grid_height = unit(5, "mm")),
               top_annotation = ha,
               border = T,
               #width = unit(100, "mm"),
               #height = unit(100, "mm"),
               use_raster = T,
               raster_quality = 5)
pdf(str_c("4.3_2_Heatmap.pdf"), width = 12, height = 6)
print(hm2)
dev.off()

# Save variants
write.table(voi.ch, file = str_c("4.3_vois_order.txt"), quote = F, col.names = F, row.names = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Visualize UMAPs and ranked VAFs #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

pdf(str_c("4.3_3_UMAPs_ranked_VAFs.pdf"))
for (v in voi.ch) {
    message(v)
    
    # Add info of variant of interest
    metadata.tib$af_voi <- af.dm[v,metadata.tib$cell]
    metadata.tib$cov_voi <- assays(maegatk.rse)[["coverage"]][as.numeric( cutf(v, d = "_") ),metadata.tib$cell]
    
    print(
        metadata.tib %>% arrange(af_voi) %>% filter(cov_voi > 5) %>%
            ggplot(aes(x = UMAP_1, y = UMAP_2, color = af_voi)) + # change to cov_voi to see coverage
            geom_point_rast(size = 1) +
            scale_color_gradientn(colors = heatcol.ch[4:10], limits = c(0,100), n.breaks = 3) +
            theme_classic() +
            theme(aspect.ratio = 1, axis.line = element_blank(), plot.title = element_text(hjust=0.5),
                  panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
            ggtitle(v)
    )
    print(
        metadata.tib %>% arrange(af_voi) %>% mutate(Rank_sorted_cells = row_number()) %>%
            ggplot(aes(x = Rank_sorted_cells, y = af_voi, color = CellType)) +
            geom_point_rast() + 
            scale_color_manual(values = mycol.ch[levels(metadata.tib$CellType)]) +
            ylab("Variant allele frequency") +
            theme(aspect.ratio = 0.5, plot.title = element_text(hjust=0.5)) +
            ggtitle(v)
    )
    metadata.tib$af_voi <- NULL
    metadata.tib$cov_voi <- NULL
}
dev.off()

# Use 683_G>A and 2593_G>A for main figure.
# To generate one figure with all remaining variants for the supplement,
# make a list of plots, then put them all on one page
p <- list()
for (v in setdiff(voi.ch, c("683_G>A", "2593_G>A"))) {
    #v <- "4104_A>G"
    metadata.tib$af_voi <- af.dm[v,metadata.tib$cell]
    metadata.tib$cov_voi <- assays(maegatk.rse)[["coverage"]][as.numeric( cutf(v, d = "_") ),metadata.tib$cell]
    
    p[[v]] <- metadata.tib %>% arrange(af_voi) %>% #filter(cov_voi > 5) %>%
        ggplot(aes(x = UMAP_1, y = UMAP_2, color = af_voi)) +
        geom_point_rast(size = 0.1, raster.dpi = 900) +
        scale_color_gradientn(colors = heatcol.ch[4:10], limits = c(0,100), n.breaks = 3) +
        theme_classic() +
        theme(aspect.ratio = 1, axis.line = element_blank(), plot.title = element_text(hjust=0.5, size = 9),
              panel.border = element_rect(colour = "black", fill=NA, size=0.5/2.13*3/2),
              legend.position = "none", axis.title = element_blank(), axis.text = element_blank(),
              axis.ticks = element_blank()) +
        ggtitle(v)
    
    metadata.tib$af_voi <- NULL
    metadata.tib$cov_voi <- NULL
}

pdf(file = str_c("4.3_4_UMAPs_combined.pdf"), height = (11-2)*1.5, width = (8.5-2)*1.5)
grid.arrange(grobs = p, ncol = 4)
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot of variant locations #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Wrangle variants
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
GenePos.tib <- tibble(Names = c("MT.ATP6", "MT.ATP8", "MT.CO1", "MT.CO2", "MT.CO3", "MT.CYB", "MT.ND1", "MT.ND2", "MT.ND3",
                                "MT.ND4", "MT.ND4L", "MT.ND5", "MT.ND6", "MT.RNR1", "MT.RNR2"),
                      start = c(8527, 8366, 5904, 7586, 9207, 14747, 3307, 4470, 10059, 10760, 10470, 12337, 14149, 648, 1671), 
                      end = c(9207, 8572, 7445, 8269, 9990, 15887, 4262, 5511, 10404, 12137, 10766, 14148, 14673, 1601, 3229))
GenePos.tib <- GenePos.tib %>% arrange(start) %>%
    mutate(mid = round((end-start)/2+start,0), ycoord = rep(c(-0.05,-0.07), length.out = 15))

pdf("4.3_5_Variant_locations.pdf", height = 3, width = 10)
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



