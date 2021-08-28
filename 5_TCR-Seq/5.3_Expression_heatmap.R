# Peter van Galen, 210628
# Evaluate subclonal structure and T-cells

# Prerequisites
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(tidyverse)
library(Seurat)
library(SummarizedExperiment)
library(readxl)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(data.table)

setwd("~/DropboxMGB/Projects/Maester/AnalysisPeter/5_TCR-Seq/")

rm(list=ls())
source("../210215_FunctionsGeneral.R")

# Functions and colors (available at https://github.com/vangalenlab/MAESTER-2021)
maester_colors.tib <- read_excel("../MAESTER_colors.xlsx")
mycol.ch <- maester_colors.tib$hex
names(mycol.ch) <- maester_colors.tib$name

# Load Seurat object containing integrated single-cell data, generated in the script 5.2_Clonal_overlap.R
seu <- readRDS("All_Seurat.rds")
seu <- ScaleData(seu)

# Load Maegtk, calculate allele frequencies (available at https://vangalenlab.bwh.harvard.edu/maester-2021/))
maegatk.rse <- readRDS("../4_CH_sample/BPDCN712_Maegatk_Final.rds")
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))*100
# Check (should be TRUE)
all(colnames(af.dm) == colnames(seu))


#~~~~~~~~~~~~~~~~~~~~~#
# Plot clones on UMAP #
#~~~~~~~~~~~~~~~~~~~~~#

# Extract data from Seurat
metadata.tib <- as_tibble(seu@meta.data, rownames = "cell")

# Choose the most interesting TCR clonotypes and MT variants empirically
metadata.tib %>% filter(TRB_CDR3 == "CASSQAGAANTEAFF") %>% dplyr::count(TNK_CellType) # Effector
metadata.tib %>% filter(TRB_CDR3 == "CASSQAGAANTEAFF") %>% dplyr::count(MT_clone) %>% arrange(desc(n)) # 1222_A>G, 2593_G>A, 779_T>C

metadata.tib %>% filter(TRB_CDR3 == "CASSLEWGNPSTYEQYF") %>%  dplyr::count(TNK_CellType) # Effector (towards Hyper)
metadata.tib %>% filter(TRB_CDR3 == "CASSLEWGNPSTYEQYF") %>% dplyr::count(MT_clone) %>% arrange(desc(n)) # 6205_G>A-9164_T>C

metadata.tib %>% filter(TRB_CDR3 == "CASSFRQGYNEQFF") %>% dplyr::count(TNK_CellType) # Hyperactivated Effector
metadata.tib %>% filter(TRB_CDR3 == "CASSFRQGYNEQFF") %>% dplyr::count(MT_clone) %>% arrange(desc(n)) # 683_G>A and 6185_T>C

metadata.tib %>% filter(TRB_CDR3 == "CASSLTGGSYNEQFF") %>% dplyr::count(TNK_CellType) # Hyperactivated Effector
metadata.tib %>% filter(TRB_CDR3 == "CASSLTGGSYNEQFF") %>% dplyr::count(MT_clone) %>% arrange(desc(n)) # 3628_A>G and 1415_G>A-9753_G>A

# Move forward with the following
tcr_clones <- c("CASSQAGAANTEAFF", "CASSLEWGNPSTYEQYF", "CASSFRQGYNEQFF", "CASSLTGGSYNEQFF")
mt_clones <- c("1222_A>G", "2593_G>A", "779_T>C", "6205_G>A-9164_T>C", "683_G>A", "6185_T>C", "3628_A>G", "1415_G>A-9753_G>A")

pdf("5.3.1_TRB-MT-highlight.pdf")

metadata.tib %>% mutate(TRB_higlight = factor(TRB_CDR3, levels = tcr_clones)) %>%
    arrange(TRB_higlight) %>%
    mutate(index = row_number()) %>% arrange(desc(index)) %>%
    ggplot(aes(x = UMAP_TNK_1, y = UMAP_TNK_2, color = TRB_higlight)) +
    geom_point(size = 0.7) +
    scale_color_manual(values = mycol.ch, na.value = "#dcdcdc") +
    theme_bw() +
    theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 3) ) )

metadata.tib %>% mutate(MT_highlight = factor(MT_clone, levels = mt_clones)) %>%
    arrange(MT_highlight) %>%
    mutate(index = row_number()) %>% arrange(desc(index)) %>%
    ggplot(aes(x = UMAP_TNK_1, y = UMAP_TNK_2, color = MT_highlight)) +
    geom_point(size = 0.7) +
    scale_color_manual(values = mycol.ch, na.value = "#dcdcdc") +
    theme_bw() +
    theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 3) ) )

dev.off()

# Subset allele frequency matrix for variants and cells of interest
af_subset.mat <- af.dm[unlist(str_split(mt_clones, "-")), filter(metadata.tib, TRB_CDR3 %in% tcr_clones, MT_clone %in% mt_clones)$cell]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Differentially expressed genes #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Determine MT clones to compare
seu.subset <- seu[,colnames(af_subset.mat)]
table(seu.subset$MT_clone) %>% sort(decreasing = T)
mt_clones_subset <- intersect(mt_clones, names(table(seu.subset$MT_clone)[table(seu.subset$MT_clone) > 10]))
seu.subset$expr_groups <- factor(seu.subset$MT_clone, levels = mt_clones_subset)
seu.subset@active.ident <- seu.subset$expr_groups
#seu.subset@active.ident <- factor(seu.subset$TRB_CDR3, levels = tcr_clones)

# Differential gene expression between mtDNA clones of interest
markerGenes <- FindAllMarkers(seu.subset, test.use = "roc", only.pos = TRUE)
markergenes.dt.ls <- lapply(split(markerGenes, f = markerGenes$cluster), function(x) data.table(x))
markergenes.dt.ls <- lapply(markergenes.dt.ls, function(x) setorder(x, -avg_logFC))
markergenes.tib <- as_tibble( do.call(cbind, lapply(markergenes.dt.ls, function(x) x$gene[1:3])) )
diff_genes <- unique(na.omit(unlist(markergenes.tib)))

# Quick visual check
DoHeatmap(seu.subset, features = diff_genes, label = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# MT, TRB, gene expression heatmap #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Customize column order
for (x in rev(strsplit(mt_clones, "-"))) {
    if (length(x) == 1) {
        af_subset.mat <- af_subset.mat[,order(-af_subset.mat[x,])]
    } else {
        af_subset.mat <- af_subset.mat[,order(-colSums(af_subset.mat[x,]))]
    }
}
# Optional: order by clonotype
#cell_order <- tibble(cell = colnames(af_subset.mat)) %>% left_join(select(metadata.tib, cell, TRB_CDR3)) %>%
#    mutate(TRB_CDR3 = factor(TRB_CDR3, levels = tcr_clones)) %>% arrange(TRB_CDR3) %>% .$cell
#af_subset.mat <- af_subset.mat[,cell_order]

# Generate gene expression matrix
plot_genes.mat <- as.matrix( GetAssayData(seu, slot = "data")[diff_genes, colnames(af_subset.mat)] )
plot_genes.mat <- plot_genes.mat - rowMeans(plot_genes.mat)
z.lim <- c(-2, 4)
plot_genes.mat[plot_genes.mat < z.lim[1]] <- z.lim[1]
plot_genes.mat[plot_genes.mat > z.lim[2]] <- z.lim[2]

# Heatmap colors
vaf.col <- colorRamp2(seq(0, 100, length.out = 9), c("#FCFCFC","#FFEDB0","#FFDF5F","#FEC510","#FA8E24","#F14C2B","#DA2828","#BE2222","#A31D1D"))
expr.col <- read_excel("~/DropboxPartners/Pipelines/AuxiliaryFiles/PopCol.xlsx", sheet = 3, col_names = F)[[1]]

# Add annotation bars
ha <- HeatmapAnnotation(vaf = t(af_subset.mat),
                        trb = seu@meta.data[colnames(af_subset.mat),"TRB_CDR3"],
                        col = list(vaf = vaf.col, trb = mycol.ch[names(mycol.ch) %in% tcr_clones]),
                        border = T)

# Plot
hm1 <- Heatmap(plot_genes.mat,
              col = colorRamp2(seq(z.lim[1], z.lim[2], length.out = 9), expr.col[3:11]),
              show_row_names = T,
              show_column_names = F,
              cluster_columns = F,
              cluster_rows = F,
              name = "expr",
              top_annotation = ha,
              border = T,
              use_raster = T)

pdf(str_c("5.3.2_Heatmap.pdf"), width = 11, height = 5)
print(hm1)
dev.off()







