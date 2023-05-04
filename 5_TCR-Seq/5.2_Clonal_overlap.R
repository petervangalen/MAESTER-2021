# Peter van Galen, 210820
# Analyze TCR-Seq data & intersect with mitochondrial variants

# Prerequisites
options(scipen = 999)

library(tidyverse)
library(Seurat)
library(readxl)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(mclust)

setwd("~/DropboxMGB/Projects/Maester/AnalysisPeter/5_TCR-Seq/")

rm(list=ls())
source("../Auxiliary_files/210215_FunctionsGeneral.R")

# Functions and colors (available at https://github.com/vangalenlab/MAESTER-2021)
maester_colors.tib <- read_excel("../MAESTER_colors.xlsx")
mycol.ch <- maester_colors.tib$hex
names(mycol.ch) <- maester_colors.tib$name

# Load data
seu.all <- readRDS("../4_CH_sample/BPDCN712_Seurat_Final.rds") # (available at https://vangalenlab.bwh.harvard.edu/maester-2021/)
seu.tnk <- readRDS("TNK_Seurat.rds") # Generated in the script 5.1_Recluster_TNK.R. Alternatively you can use this (available at https://vangalenlab.bwh.harvard.edu/maester-2021/): seu.tnk <- readRDS("BPDCN712_Seurat_with_TCR.rds")

# Import TCR data. Note that there are no duplicate cell barcodes because Duncan processed the data.
tcr.tib <- read_csv("Duncan_210805/20210805_maester2_vfilter9.csv")[,-1]
# Filter for high-quality cells
tcr.tib <- tcr.tib %>% mutate(cell = str_c(BC, "-1")) %>% select(-BC) %>% filter(cell %in% colnames(seu.all))

# Some numbers for the Methods text
T.CTL <- colnames(seu.tnk)[grepl("T", seu.tnk$CellType)]
length(T.CTL)
sum(unique(filter(tcr.tib, grepl("TRB", TCR_Recovery))$cell) %in% T.CTL)
sum(unique(filter(tcr.tib, grepl("TRA", TCR_Recovery))$cell) %in% T.CTL)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Determine & plot top TRB clones #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Only keep cells with TRB calls found in at least 10 cells (this is very stringent)
tcr_filter.tib <- tcr.tib %>% filter(!is.na(TRB_CDR3)) %>%
    group_by(TRB_CDR3) %>% mutate(n_b = n()) %>% filter(n_b > 10) %>% ungroup() %>%
    arrange(desc(n_b), TRB_CDR3, TRA_CDR3) %>% 
    mutate(TRB_CDR3 = factor(TRB_CDR3, levels = unique(TRB_CDR3))) %>%
    select(cell, TRB_CDR3, TRA_CDR3)

# Only keep TRA clones called in at least 10 cells
top_a <- tcr_filter.tib %>% dplyr::count(TRA_CDR3) %>% filter(n > 10) %>%
    arrange(desc(n)) %>% na.omit %>% .$TRA_CDR3
tcr_filter.tib <- tcr_filter.tib %>% mutate(TRA_CDR3 = factor(TRA_CDR3, levels = top_a)) %>%
    arrange(TRB_CDR3, TRA_CDR3)

# Combine with other metadata (mainly for UMAP coordinates)
metadata.tib <- as_tibble(seu.all@meta.data, rownames = "cell")
metadata.tib <- left_join(metadata.tib, tcr_filter.tib, by = "cell")
metadata.tib <- left_join(metadata.tib, select(as_tibble(seu.tnk@meta.data, rownames = "cell"),
    cell, UMAP_TNK_1, UMAP_TNK_2, TNK_CellType), by = "cell")

# Visualize
pdf("5.2.1_TRB_clone_UMAPs.pdf")

metadata.tib %>% filter(! is.na(TNK_CellType)) %>%
    ggplot(aes(x = UMAP_TNK_1, y = UMAP_TNK_2, color = TNK_CellType)) +
    geom_point(size = 0.7) +
    scale_color_manual(values = mycol.ch[levels(metadata.tib$TNK_CellType)]) +
    ggtitle("TNK clusters") +
    theme_bw() +
    theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 3) ) )

for (z in levels(metadata.tib$TRB_CDR3)) {
    print(
        metadata.tib %>% filter(! is.na(TNK_CellType)) %>%
            mutate(current_clone = replace_na(TRB_CDR3 == z, F)) %>%
            arrange(current_clone) %>%
            ggplot(aes(x = UMAP_TNK_1, y = UMAP_TNK_2, color = current_clone)) +
            geom_point(size = 0.7) +
            scale_color_manual(values = c("lightgrey", "purple")) +
            ggtitle(label = z) +
            theme_bw() +
            theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            guides(color = guide_legend(override.aes = list(size = 3) ) )
    )
}

dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Overlap between TRA and TRB #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Wrangle
heatmap.tib <- metadata.tib %>% group_by(TRB_CDR3, TRA_CDR3) %>% mutate(clone_size = n()) %>% #filter(clone_size > 10) %>%
    select(TRB_CDR3, TRA_CDR3, CellType) %>% na.omit() %>% arrange(TRB_CDR3, TRA_CDR3, CellType)
heatmap.tib <- heatmap.tib %>%
    mutate(TRB_CDR3 = factor(TRB_CDR3, levels = unique(TRB_CDR3))) %>%
    mutate(TRA_CDR3 = factor(TRA_CDR3, levels = unique(TRA_CDR3)))

# Assess agreement
adjustedRandIndex(heatmap.tib$TRB_CDR3, heatmap.tib$TRA_CDR3) # 0.8411752

# Plot
ht <- Heatmap(t(as.matrix(heatmap.tib)),
              col = c(mycol.ch),
              na_col = "#FFFFFF",
              cluster_rows = F,
              cluster_columns = F,
              show_row_names = T,
              show_column_names = T,
              show_heatmap_legend = F,
              column_title = str_c(nrow(heatmap.tib), " cells"))

# Assemble three legends
TRB.lgd <- Legend(at = as.character(unique(heatmap.tib$TRB_CDR3)),
                  legend_gp = gpar(fill = mycol.ch[as.character(unique(heatmap.tib$TRB_CDR3))]), title = "TRB variable region")
TRA.lgd <- Legend(at = as.character(unique(heatmap.tib$TRA_CDR3)),
                  legend_gp = gpar(fill = mycol.ch[as.character(unique(heatmap.tib$TRA_CDR3))]), title = "TRA variable region")
CellType.lgd <- Legend(at = names(mycol.ch[1:14]), legend_gp = gpar(fill = mycol.ch[1:14]), title = "Cell type")
legend_pack <- packLegend(TRB.lgd, TRA.lgd, CellType.lgd, direction = "horizontal" )

# Plot
grobs.ls <- list(grid.grabExpr(draw(ht)), grid.grabExpr(draw(legend_pack)))

pdf("5.2.2_TRB-TRA.pdf", width = 8, height = 9)
grid.arrange(grobs = grobs.ls, ncol = 1, nrow = 2, heights = c(1,3))
dev.off()


#~~~~~~~~~~~~~~~~#
# Load MT clones #
#~~~~~~~~~~~~~~~~#

# Import clonal data (available at https://github.com/vangalenlab/MAESTER-2021)
positive_cells.tib <- read_tsv("../4_CH_sample/4.4_positive_cells.txt")
positive_cells.tib <- positive_cells.tib %>% mutate(clone = factor(clone, levels = unique(clone)),
                                                    variant = factor(variant, levels = unique(variant)),
                                                    CellType = factor(CellType, levels = maester_colors.tib$name[1:14]))

# A small number of cells is assigned to multiple MT clones
positive_cells.tib %>% group_by(cell) %>%
    group_modify(~ {.x %>% mutate(n_variants = length(unique(variant)))}) %>%
    group_modify(~ {.x %>% mutate(n_clones = length(unique(clone)))}) %>%
    filter(n_clones > 1) %>% arrange(desc(n_clones), desc(n_variants), cell, desc(cov))
# Only the MT clone with the highest coverage is kept to facilitate downstream processing.
filtered_cells.tib <- positive_cells.tib %>% mutate(index = row_number()) %>% group_by(cell) %>% 
    group_modify(~ {.x %>% arrange(desc(cov))}) %>%
    filter(!duplicated(cell)) %>%
    arrange(clone, index) %>%
    dplyr::rename(MT_clone = clone)

# Combine with Seurat metadata
metadata.tib <- left_join(metadata.tib, select(filtered_cells.tib, cell, MT_clone), by = "cell")

pdf("5.2.3_MT_clone_UMAPs.pdf")

metadata.tib %>% filter(! is.na(TNK_CellType)) %>%
    ggplot(aes(x = UMAP_TNK_1, y = UMAP_TNK_2, color = TNK_CellType)) +
    geom_point(size = 0.7) +
    scale_color_manual(values = mycol.ch[levels(metadata.tib$TNK_CellType)]) +
    ggtitle("TNK clusters") +
    theme_bw() +
    theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 3) ) )

for (z in levels(metadata.tib$MT_clone)) {
    print(
        metadata.tib %>% filter(! is.na(TNK_CellType)) %>%
            mutate(current_clone = replace_na(MT_clone == z, F)) %>%
            arrange(current_clone) %>%
            ggplot(aes(x = UMAP_TNK_1, y = UMAP_TNK_2, color = current_clone)) +
            geom_point(size = 0.7) +
            scale_color_manual(values = c("lightgrey", "purple")) +
            ggtitle(label = z) +
            theme_bw() +
            theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            guides(color = guide_legend(override.aes = list(size = 3) ) )
    )
}
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Evaluate TRB / MT overlap #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Determine largest MT clone in TRB clones and vice versa, and the size of clones defined by both
MT_TRB.tib <- metadata.tib %>% select(MT_clone, TRB_CDR3, CellType) %>%
    group_by(TRB_CDR3) %>%
    group_modify(~ {.x %>% mutate(largest_mt = max(table(MT_clone)))}) %>%
    ungroup %>% group_by(MT_clone) %>%
    group_modify(~ {.x %>% mutate(largest_trb = max(table(TRB_CDR3)))}) %>%
    ungroup %>% group_by(TRB_CDR3, MT_clone) %>%
    mutate(n_both = n()) %>%
    ungroup

# Select only TRB clones that have an MT clone of >5 cells and vice versa
MT_TRB_filter.tib <- MT_TRB.tib %>% filter(largest_mt > 5, largest_trb > 5)

# Reorder TCR clones manually (based on how it looks in the confusion matrix later)
TRB_clone_order <- c("CASSLVEEKLFF", "CASSFRQGYNEQFF", "CASSLEWGNPSTYEQYF", "CASSLTGGSYNEQFF", "CASSRYRGGTEAFF", "CASSPFEETQYF",
                     "CSVEDLGMDSPSYNEQFF", "CATSGGWGGETQYF", "CASSQAGAANTEAFF", "CASSQVGHSADTQYF", "CASSLTSGDPTDTQYF")
MT_TRB_filter.tib %>% group_by(TRB_CDR3) %>% dplyr::count() %>% mutate(shown = TRB_CDR3 %in% TRB_clone_order)
MT_TRB_filter.tib <- MT_TRB_filter.tib %>% mutate(TRB_CDR3 = factor(TRB_CDR3, levels = TRB_clone_order))

# Wrangle
heatmap.tib <- MT_TRB_filter.tib %>% select(MT_clone, TRB_CDR3) %>% na.omit %>% arrange(TRB_CDR3, MT_clone)
heatmap.tib <- heatmap.tib %>%
    mutate(MT_clone = factor(MT_clone, levels = intersect(levels(MT_clone), heatmap.tib$MT_clone))) %>%
    mutate(TRB_CDR3 = factor(TRB_CDR3, levels = intersect(levels(TRB_CDR3), heatmap.tib$TRB_CDR3)))

# Assess agreement
adjustedRandIndex(heatmap.tib$MT_clone, heatmap.tib$TRB_CDR3) # .7394777
write_tsv(heatmap.tib, file = "5.2_Clonal_overlap.txt")

# Plot
ht <- Heatmap(t(as.matrix(heatmap.tib)),
              col = mycol.ch,
              na_col = "#FFFFFF",
              cluster_rows = F,
              cluster_columns = F,
              show_row_names = T,
              show_column_names = T,
              show_heatmap_legend = F,
              column_title = str_c(nrow(heatmap.tib), " cells")
              )

# Assemble legends
TRB.lgd <- Legend(at = as.character(unique(heatmap.tib$TRB_CDR3)),
                  legend_gp = gpar(fill = mycol.ch[as.character(unique(heatmap.tib$TRB_CDR3))]), title = "TRB variable region")
MT.lgd <- Legend(at = as.character(unique(heatmap.tib$MT_clone)),
                 legend_gp = gpar(fill = mycol.ch[as.character(unique(heatmap.tib$MT_clone))]), title = "MT clone")
legend_pack <- packLegend(MT.lgd, TRB.lgd, direction = "horizontal" )
# Ensure there are no duplicate colors in the legend
duplicated(c(TRB.lgd@grob$children[[2]]$children[[2]]$gp$fill,
             MT.lgd@grob$children[[2]]$children[[2]]$gp$fill))

# Plot
grobs.ls <- list(grid.grabExpr(draw(ht)), grid.grabExpr(draw(legend_pack)))

pdf("5.2.4_MT-TRB.pdf", width = 8, height = 6)
grid.arrange(grobs = grobs.ls, ncol = 1, nrow = 2, heights = c(1,3))
dev.off()


#~~~~~~~~~~~~~~~~~~#
# Confusion matrix #
#~~~~~~~~~~~~~~~~~~#

# Normalize clone sizes to cell number with the mitochondrial mutation
plot.tib <- MT_TRB_filter.tib %>% select(MT_clone, TRB_CDR3, n_both) %>% na.omit %>%
    group_by(MT_clone) %>% mutate(confusion_color = n_both / n())

pdf("5.2.5_ConfusionPlot_tidy.pdf")

plot.tib %>% mutate(TRB_CDR3 = factor(TRB_CDR3, levels = rev(levels(TRB_CDR3)))) %>%
    ggplot(aes(x = MT_clone, y = TRB_CDR3, fill = confusion_color)) +
    geom_tile() +
    scale_fill_gradientn(colors = c("#ffffff", "#ff0000", "#dc143c")) +
    #geom_jitter(width = 0.4, height = 0.4, size = 0.1) +
    coord_fixed() +
    geom_hline(yintercept = seq(0.5, 40, 1), size = 0.2) +
    geom_vline(xintercept = seq(0.5, 40, 1), size = 0.2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          legend.background = element_blank(),
          axis.ticks = element_blank())

dev.off()


# If you thought that was confusing, here is a base R approach
# Generate a matrix with group sizes
confusion.mat <- matrix(NA, nrow = length(levels(metadata.tib$TRB_CDR3)), ncol = length(levels(metadata.tib$MT_clone)))
rownames(confusion.mat) <- levels(metadata.tib$TRB_CDR3)
colnames(confusion.mat) <- levels(metadata.tib$MT_clone)
for (a in rownames(confusion.mat)) {
    for (b in colnames(confusion.mat)) {
        confusion.mat[a,b] <- metadata.tib %>% filter(TRB_CDR3 == a, MT_clone == b) %>% nrow
    }
}

# Select rows and columns with at least 5 in one column or one row
confusion_subset.mat <- confusion.mat[rowSums(confusion.mat > 5) > 0,
                                      colSums(confusion.mat > 5) > 0]
confusion.mat %>% dim
confusion_subset.mat %>% dim

# Normalize to column sums
confusion_norm.mat <- sweep(confusion_subset.mat, 2, colSums(confusion_subset.mat), "/")*100
confusion_norm.mat <- confusion_norm.mat[TRB_clone_order,]

# Plot
pdf("5.2.5_ConfusionPlot_base.pdf")

Heatmap(confusion_norm.mat,
        col = colorRamp2(breaks = c(0, 50, 100), colors = c("#ffffff", "#ff0000", "#dc143c")),
        cluster_rows = F,
        cluster_columns = F,
        rect_gp = gpar(col = "black", lwd = 1),
        row_title = "TCR beta chain",
        row_names_side = "left",
        column_title = "MT clone",
        column_title_side = "bottom",
        height = unit(nrow(confusion_norm.mat)*0.6, "cm"),
        width = unit(ncol(confusion_norm.mat)*0.6, "cm")
        )

dev.off()


# Save for the next script
all( seu.all@meta.data %>% rownames == metadata.tib$cell )
setdiff(colnames(metadata.tib), colnames(seu.all@meta.data))
seu.all@meta.data$UMAP_TNK_1 <- metadata.tib$UMAP_TNK_1
seu.all@meta.data$UMAP_TNK_2 <- metadata.tib$UMAP_TNK_2
seu.all@meta.data$TNK_CellType <- metadata.tib$TNK_CellType
seu.all@meta.data$TRB_CDR3 <- metadata.tib$TRB_CDR3
seu.all@meta.data$TRA_CDR3 <- metadata.tib$TRA_CDR3
seu.all@meta.data$MT_clone <- metadata.tib$MT_clone

saveRDS(seu.all, file = "BPDCN712_Seurat_with_TCR.rds")



