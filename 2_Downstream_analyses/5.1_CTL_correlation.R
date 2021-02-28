# Peter van Galen, 210201
# Correlation of cells within (intraclonal) and between (interclonal) CTL clones

# Prerequisites
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(tidyverse)
library(Seurat)
library(readxl)
library(ComplexHeatmap)
library(circlize) # for colorRamp2

rm(list=ls())
setwd("~/DropboxPartners/Projects/Maester/AnalysisPeter/5_CTL_correlation")

# Functions and colors (available at https://github.com/vangalenlab/MAESTER-2021)
source("../210215_FunctionsGeneral.R")
popcol.tib <- read_excel("../MAESTER_colors.xlsx")
mycol.ch <- popcol.tib$hex
names(mycol.ch) <- popcol.tib$name

# Load Seurat object (available at https://vangalenlab.bwh.harvard.edu/maester-2021/)
seu <- readRDS("../4_CH_sample/BPDCN712_Seurat.rds")

# Subset for CTLs & identify variable genes
CTL.seu <- subset(seu, subset = CellType == "CTL")
CTL.seu <- FindVariableFeatures(CTL.seu, nfeatures = 100)
tibble(Variable_features = VariableFeatures(CTL.seu)) %>% write_tsv("VariableFeatures.txt")

# Cell IDs
positive_cells.tib <- read_tsv(file = "../4_CH_sample/210204_positive_cells.txt")
positive_cells.tib <- positive_cells.tib %>% filter(!duplicated(cell)) %>% mutate(clone = factor(clone, levels = unique(clone)),
                                                    variant = factor(variant, levels = unique(variant)),
                                                    CellType = factor(CellType, levels = popcol.tib$name[1:14]))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Heatmap of 1,044 x 1,044 CTLs ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Make a tibble with CTLs and clones
clones.sorted <- positive_cells.tib %>% group_by(clone) %>% summarize(n = n()) %>% arrange(desc(n)) %>% .$clone
CTL.tib <- filter(positive_cells.tib, CellType == "CTL")

# Get expression data of cells in clones
expr.mat <- as.matrix( GetAssayData(CTL.seu, slot = "data") )
expr.mat <- expr.mat[VariableFeatures(CTL.seu),CTL.tib$cell]

ha <- HeatmapAnnotation(df = data.frame(dplyr::select(CTL.tib, clone)), col = list(clone = mycol.ch))

cor.mat <- cor(expr.mat)

pdf("1_Heatmap.pdf")
Heatmap(cor.mat,
        col = colorRamp2(c(0,0.6,1), c("blue", "white", "red")),
        cluster_rows = F,
        cluster_columns = F,
        show_row_names = F,
        show_column_names = F,
        top_annotation = ha,
        width = unit(100, "mm"),
        height = unit(100, "mm"),
        use_raster = T)
dev.off()


#~~~~~~~~~~~~~~~~#
#### Averages ####
#~~~~~~~~~~~~~~~~#

# Prepare columns to fill in
CTL_clones.tib <- CTL.tib %>% group_by(clone) %>% summarize(n = n()) %>%
    mutate(intra_cor = as.numeric(NA), inter_cor = as.numeric(NA))

# Expression matrix
expr.mat <- as.matrix(GetAssayData(CTL.seu, slot = "data"))[VariableFeatures(CTL.seu),]

for (i in CTL_clones.tib$clone) {
    #i <- "4104_A>G"
    message(i)
    current.cells <- CTL.tib %>% filter(clone == i) %>% .$cell
    other.cells <- setdiff(colnames(CTL.seu), current.cells)
    #other.cells <- CTL.tib %>% filter(clone != i) %>% .$cell
    
    # Calculate gene expression correlation of every intraclonal CTL with every other CTL of the same clone
    cor_intra.ls <- lapply(current.cells, function(x) { cor(expr.mat[,x], expr.mat[,setdiff(current.cells, x)]) } )
    names(cor_intra.ls) <- current.cells

    # Calculate gene expression correlation of every CTL with every CTL in other clones
    cor_inter.ls <- lapply(current.cells, function(x) { cor(expr.mat[,x], expr.mat[,other.cells]) } )
    names(cor_inter.ls) <- current.cells

    CTL_clones.tib[which(CTL_clones.tib$clone == i),"intra_cor"] <- mean(unlist(cor_intra.ls))
    CTL_clones.tib[which(CTL_clones.tib$clone == i),"inter_cor"] <- mean(unlist(cor_inter.ls))
}

pdf("2_Correlation_dotplot.pdf", height = 7, width = 5)

CTL_clones.tib %>% filter(n > 10) %>%
    pivot_longer(cols = c("intra_cor", "inter_cor"), values_to = "Correlation") %>%
    mutate(name = gsub("intra_cor", "intraclonal", gsub("inter_cor", "interclonal", name))) %>%
    mutate(name = factor(name, levels = c("intraclonal", "interclonal"))) %>%
    ggplot(aes(x = name, y = Correlation, color = clone, group = clone)) +
    geom_path(size = 1) +
    geom_point(size = 4) +
    scale_color_manual(values = mycol.ch) +
    xlab("") +
    theme_classic() +
    labs(size = 10) +
    theme(aspect.ratio = 3, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 14),
          axis.title.y = element_text(color = "black", size = 14),
          axis.text.y = element_text(color = "black", size = 14),
          axis.ticks = element_line(color = "black"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10))

dev.off()

# Significant?
t.test(filter(CTL_clones.tib, n >10)$intra_cor, filter(CTL_clones.tib, n > 10)$inter_cor, paired = T, var.equal = T)$p.value

