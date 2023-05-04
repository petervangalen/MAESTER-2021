# Peter van Galen, 210817
# Identify informative variants in glioblastoma sample

#~~~~~~~~~~~~~~~#
# Prerequisites #
#~~~~~~~~~~~~~~~#

options(scipen = 999)

library(tidyverse)
library(data.table)
library(SummarizedExperiment)
library(readxl)
library(ggrastr)
library(Matrix)
library(GGally) # for ggpairs
library(ggrepel)
library(gridExtra)

rm(list=ls())
setwd("~/DropboxMGB/Projects/Maester/AnalysisPeter/8_Glioblastoma")

## Functions & colors (available at https://github.com/vangalenlab/MAESTER-2021)
source("../Auxiliary_files/210215_FunctionsGeneral.R")
heat.col <- read_excel("../MAESTER_colors.xlsx", sheet = 2, col_names = "heatcol")$heatcol
celltype.col <- c(AC = "#ffebcd", MES1 = "#ffd700", MES2 = "#bc8f8f", NPC1 = "#daa520", NPC2 = "#ff4500",
                  OPC = "#7fff00", Oligo = "#4682b4", Tumor = "#00ff7f", Malignant = "#b0e0e6", Cycling = "#5f9ea0",
                  Stromal = "#ff00ff", Myeloid = "#a9a9a9", Neutrophils = "#66cdaa", Mast = "#000000", pDCs = "#00bfff",
                  Tcells = "#696969", Bcells = "#008080", Plasma = "#9932cc")
subset.col <- c(A = "#ff69b4", C = "#4b0082", PBMC = "#66cdaa", NonMalignant = "#d3d3d3")


#~~~~~~~~~~~#
# Load data #
#~~~~~~~~~~~#

start_time <- Sys.time() # This section takes about 5 minutes

# Load locally saved data. This data is available at https://vangalenlab.bwh.harvard.edu/maester-2021/
maegatk_A.full <- readRDS(file = "../1_MT_Coverage/SW_MGH252_A_mr3_maegatk.rds")
maegatk_C.full <- readRDS(file = "../1_MT_Coverage/SW_MGH252_C_mr3_maegatk.rds")
maegatk_PBMC.full <- readRDS(file = "../1_MT_Coverage/SW_MGH252_PBMC_mr3_maegatk.rds")

# Wrangle & make cell barcodes match. This data is available at https://github.com/vangalenlab/MAESTER-2021/
metadata.full <- read_tsv("Chadi email 210721/MGH252_A_C_PBMC_MetaData.txt")
umap.full <- read_tsv("Chadi email 210721/MGH252_A_C_PBMC_UMAP.txt") 
metadata.full <- metadata.full %>% left_join(umap.full) %>%
    mutate(Site = gsub("[0-9]", "", cutf(...1, d = "\\.|_", f = 4)), .before = 2) %>%
    mutate(cell = str_c(cutf(...1, d = "\\.", f = 2), "-", Site), .before = 2) %>%
    select(cell, Site, CellType, UMAP_1, UMAP_2)
colnames(maegatk_A.full) <- gsub("-1", "-A", colnames(maegatk_A.full))
colnames(maegatk_C.full) <- gsub("-1", "-C", colnames(maegatk_C.full))
colnames(maegatk_PBMC.full) <- gsub("-1", "-PBMC", colnames(maegatk_PBMC.full))

# Intersect barcodes for maegatk and metadata
metadata.full <- metadata.full %>% group_by(cell) %>% filter(n() == 1) %>% ungroup
cells.ls <- list(A = intersect(metadata.full$cell, colnames(maegatk_A.full)),
                 C = intersect(metadata.full$cell, colnames(maegatk_C.full)),
                 PBMC = intersect(metadata.full$cell, colnames(maegatk_PBMC.full)))
metadata.tib <- metadata.full %>% filter(cell %in% unlist(cells.ls))
maegatk <- cbind(maegatk_A.full[,cells.ls$A], maegatk_C.full[,cells.ls$C], maegatk_PBMC.full[,cells.ls$PBMC])
# Clear some space
rm(list=setdiff(ls()[grepl("maegatk", ls())], "maegatk"))

# Wrangle more
metadata.tib <- metadata.tib %>% mutate(CellType = factor(CellType, levels = c("AC", "MES1", "MES2", "NPC1", "NPC2", "OPC", "Oligo",
                                                                               "Tumor", "Malignant", "Cycling", "Stromal", "Myeloid", 
                                                                               "Neutrophils", "Mast", "pDCs", "Tcells", "Bcells", "Plasma")))

# Visualize cell types
pdf("8.1_1_UMAP.pdf")
metadata.tib %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = CellType)) +
    geom_point_rast(size = 0.2) +
    scale_color_manual(values = celltype.col) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    theme_bw() +
    theme(aspect.ratio = 1, panel.grid = element_blank())
dev.off()

# Calculate allele frequencies
af.dm <- data.matrix(computeAFMutMatrix(maegatk))*100
all(colnames(af.dm) == metadata.tib$cell)

Sys.time() - start_time

write_tsv(metadata.tib, file = "MGH252_metadata.txt")
saveRDS(maegatk, file = "MGH252_maegatk.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Which cell types to compare #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

malignant.ch <- c("AC", "MES1", "MES2", "NPC1", "NPC2", "OPC", "Tumor", "Malignant", "Cycling")
nonmalignant.ch <- c("Oligo", "Stromal", "Myeloid",  "Neutrophils", "Mast", "pDCs", "Tcells", "Bcells", "Plasma")
all(metadata.tib$CellType %in% c(malignant.ch, nonmalignant.ch))

CellSubsets.ls <- list(unionCells = metadata.tib$cell,
                       NonMalignant = filter(metadata.tib, Site != "PBMC", CellType %in% nonmalignant.ch)$cell,
                       A = filter(metadata.tib, Site == "A", CellType %in% malignant.ch)$cell,
                       C = filter(metadata.tib, Site == "C", CellType %in% malignant.ch)$cell,
                       PBMC = filter(metadata.tib, Site == "PBMC")$cell)
lengths(CellSubsets.ls)

# Visualize
CellSubsets.tib <- tibble(cell = unlist(CellSubsets.ls), CellSubset = gsub("[0-9]", "", names(unlist(CellSubsets.ls))))

pdf(str_c("8.1_2_CellSubsets.pdf"))
metadata.tib %>% left_join(filter(CellSubsets.tib, CellSubset != "unionCells")) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = CellSubset)) +
    geom_point_rast(size = 0.2) +
    scale_color_manual(values = subset.col) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    theme_bw() +
    theme(aspect.ratio = 1, panel.grid = element_blank())
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Collect variant information #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

start_time <- Sys.time() # This section takes about 8 minutes

# Get the mean coverage for every cell subset
mean_af.ls <- lapply(CellSubsets.ls, function(x) rowMeans(af.dm[,x]))
mean_cov.ls <- lapply(CellSubsets.ls, function(x) rowMeans(assays(maegatk)[["coverage"]][,x])[as.numeric(cutf(rownames(af.dm), d = "_"))])
names(mean_af.ls) <- paste0("mean_af.", names(mean_af.ls))
names(mean_cov.ls) <- paste0("mean_cov.", names(mean_cov.ls))

# Get the quantiles of the VAFs of each variant in each cell subset
quantiles <- c("q99" = 0.99, "q999" = 0.999)
quantiles.ls <- lapply(quantiles, function(x) lapply(CellSubsets.ls, function(y) apply(af.dm[,y], 1, quantile, x) ))

# Get the mean quality for each variant. This is 100x faster on regular matrices than the dgCMatrix.
assays.ls <- lapply(maegatk@assays@data, function(x) as.matrix(x)) # Not sure why it's not $data as for the CH sample
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

# Collect all information in a tibble
vars.tib <- as_tibble(do.call(cbind, c(mean_af.ls, mean_cov.ls, unlist(quantiles.ls, recursive = F))), rownames = "var")
vars.tib <- add_column(vars.tib, quality = qual.num, .before = 2)

Sys.time() - start_time

# Save for fast loading next time
write_tsv(vars.tib, "8.1_MGH252_variants.txt")
vars.tib <- read_tsv("8.1_MGH252_variants.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Select variants of interest #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# This yields 19 variants of which 89% transitions
voi_A.ch <- vars.tib %>% filter(mean_cov.unionCells > 5, quality >= 30) %>%
    filter(q99.A > 5, q99.A > 10*q99.C, q99.A > 10*q99.PBMC) %>% .$var
voi_C.ch <- vars.tib %>% filter(mean_cov.unionCells > 5, quality >= 30) %>%
    filter(q99.C > 5, q99.C > 10*q99.A, q99.C >10*q99.PBMC) %>% .$var
voi_PBMC.ch <- vars.tib %>% filter(mean_cov.unionCells > 5, quality >= 30) %>%
    filter(q99.PBMC > 5, q99.PBMC > 10*q99.A, q99.PBMC > 10*q99.C) %>% .$var
voi.ch <- c(voi_A.ch, voi_C.ch, voi_PBMC.ch)

# Assess transitions vs. transversions
str_view(voi.ch, "G>A|A>G|C>T|T>C"); length(voi.ch); mean( str_count(voi.ch, "G>A|A>G|C>T|T>C") )

# Visualize pairs. Use png instead of pdf b/c there are too many data points
png("8.1.3_Scatter_pairs.png", width = 1500, height = 1400, res = 300, pointsize = 0.6)
print(
vars.tib %>% filter(mean_cov.unionCells > 5, quality >= 30) %>%
    select(var, contains("q99.")) %>%
    mutate(cell_subset = case_when(
        var %in% voi_A.ch ~ "A",
        var %in% voi_C.ch ~ "C",
        var %in% voi_PBMC.ch ~ "PBMC")) %>%
    arrange(!is.na(cell_subset), cell_subset) %>%
    rename(q99.A = "A", q99.C = "C", q99.PBMC = "PBMC") %>%
    ggpairs(lower = list(continuous = wrap("points", size = 0.8)),
        diag = list(continuous = wrap("barDiag", binwidth = 10, fill = "grey")),
        upper = list(continuous = wrap("points", size = 0.8)),
        columns = c("A", "C", "PBMC"),
        aes(color = cell_subset)) +
    scale_color_manual(values = subset.col, na.value = "#e6e6fa") +
    #geom_text_repel(mapping = aes())
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black"),
          axis.ticks = element_line(color = "black"),
          panel.border = element_rect(colour = "black"))
)
dev.off()


# Or a simpler visualization:
gg.list <- list()

for (l in list(c("A", "C"), c("A", "PBMC"), c("C", "PBMC"))) {
    #l <- c("A", "C")
    vars2.tib <- vars.tib %>% filter(mean_cov.unionCells > 5, quality >= 30) %>%
        select(var, contains("q99.")) %>%
        mutate(cell_subset = case_when(
            var %in% voi_A.ch ~ "A",
            var %in% voi_C.ch ~ "C",
            var %in% voi_PBMC.ch ~ "PBMC")) %>%
        arrange(!is.na(cell_subset), cell_subset) %>%
        rename(q99.A = "A", q99.C = "C", q99.PBMC = "PBMC")

    g <- vars2.tib %>%
        ggplot(aes_string(x = l[1], y = l[2], color = "cell_subset")) +
        geom_point_rast() +
        geom_text_repel(data = filter(vars2.tib, var %in% c("2817_G>A", "15635_T>C", "2056_G>A")),
                        aes_string(x = l[1], y = l[2], label = "var"), show.legend = FALSE) +
        scale_color_manual(values = subset.col, na.value = "#c0c0c0") +
        theme_bw() +
        theme(aspect.ratio = 1,
              panel.grid = element_blank(),
              axis.text = element_text(color = "black"),
              axis.ticks = element_line(color = "black"),
              panel.border = element_rect(colour = "black"))
    
    gg.list[[length(gg.list) + 1]] <- g
    }

pdf("8.1_3_Scatter_simple.pdf", width = 10)
grid.arrange(gg.list[[1]], gg.list[[2]], gg.list[[3]], nrow = 1)
dev.off()


#~~~~~~~~~~~~~~~~~~#
#### Plot UMAPs ####
#~~~~~~~~~~~~~~~~~~#

pdf(str_c("8.1_4_MGH252_UMAPs.pdf"))

for (v in voi.ch) {
    # v <- voi.ch[1]
    message(v)
    
    # Add info of variant of interest
    metadata.tib$af_voi <- af.dm[v,metadata.tib$cell]
    metadata.tib$cov_voi <- assays(maegatk)[["coverage"]][as.numeric( cutf(v, d = "_") ),metadata.tib$cell]
    
    print(
        metadata.tib %>% arrange(af_voi) %>% filter(cov_voi >= 5) %>%
            ggplot(aes(x = UMAP_1, y = UMAP_2, color = af_voi)) +
            geom_point_rast(size = 1) +
            scale_color_gradientn(colors = heat.col[4:10], n.breaks = 3) +
            theme_classic() +
            theme(aspect.ratio = 1, axis.line = element_blank(), plot.title = element_text(hjust=0.5),
                  panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
            ggtitle(v)
    )
    print(
        metadata.tib %>% arrange(af_voi) %>% mutate(Rank_sorted_cells = row_number()) %>%
            ggplot(aes(x = Rank_sorted_cells, y = af_voi, color = CellType, shape = cov_voi >= 5)) +
            scale_color_manual(values = celltype.col) +
            scale_shape_manual(values=c(4, 16))+
            geom_point_rast() +
            ylab("Variant allele frequency") +
            theme(aspect.ratio = 0.5, plot.title = element_text(hjust=0.5)) +
            ggtitle(v)
    )
    metadata.tib$af_voi <- NULL
    metadata.tib$cov_voi <- NULL
}

dev.off()

# Save variants
write_tsv(tibble(var = voi.ch), file = str_c("8.1_vois.txt"))


### Stats for legend

# Make a list of cell ids that are positive for each of the variants of interest
positive_cells.ls <- list()

for (v in voi.ch) {
    message(v)
    # Determine cells with an appreciable VAF
    current_cells.ch <- colnames(af.dm)[af.dm[v,]>1]
    # Save cell IDs for positive cells
    positive_cells.ls[[v]] <- current_cells.ch
}

tracked_cells <- unlist(positive_cells.ls) %>% unique %>% length
total_cells <- ncol(af.dm)
tracked_cells / total_cells * 100


