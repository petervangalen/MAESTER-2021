# Peter van Galen, 211111
# Test different parameters to select informative variants from the clonal hematopoiesis sample.


#~~~~~~~~~~~~~~~~~~~~~#
#### Prerequisites ####
#~~~~~~~~~~~~~~~~~~~~~#

options(scipen = 999)

library(tidyverse)
library(SummarizedExperiment)
library(Seurat)
library(Matrix)
library(readxl)
library(ComplexHeatmap)
library(circlize) # For heatmap colors
library(magick) # For heatmap rasterization
library(mclust) # For adjusted rand index

rm(list=ls())
setwd("~/DropboxMGB/Projects/Maester/AnalysisPeter/4_CH_sample")

# Function & colors
source("../210215_FunctionsGeneral.R")
popcol.tib <- read_excel("../MAESTER_colors.xlsx")
mycol.ch <- popcol.tib$hex
names(mycol.ch) <- popcol.tib$name

# Load Seurat object (available at https://vangalenlab.bwh.harvard.edu/maester-2021/)
seu <- readRDS("../5_TCR-Seq/BPDCN712_Seurat_with_TCR.rds")
metadata.tib <- as_tibble(seu@meta.data, rownames = "cell")

# Load Maegtk, calculate allele frequencies (https://vangalenlab.bwh.harvard.edu/maester-2021/))
maegatk.rse <- readRDS("../4_CH_sample/BPDCN712_Maegatk_Final.rds")
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))*100


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Collect information for each variant #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Get the mean quality for each variant (2.5 minutes)
start_time <- Sys.time()
assays.ls <- lapply(maegatk.rse@assays$data, function(x) as.matrix(x))
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

# Make tibble
vars.tib <- tibble(var = rownames(af.dm),
                   mean_af = rowMeans(af.dm),
                   mean_cov = rowMeans(assays(maegatk.rse)[["coverage"]])[as.numeric(cutf(rownames(af.dm), d = "_"))],
                   quality = qual.num)

# Calculate the number of cells that exceed VAF thresholds 0, 1, 5, 10, 50 (3 minutes)
start_time <- Sys.time()
vars.tib <- vars.tib %>%
    mutate(n0 = apply(af.dm, 1, function(x) sum(x == 0))) %>%
    mutate(n1 = apply(af.dm, 1, function(x) sum(x > 1))) %>%
    mutate(n5 = apply(af.dm, 1, function(x) sum(x > 5))) %>%
    mutate(n10 = apply(af.dm, 1, function(x) sum(x > 10))) %>%
    mutate(n50 = apply(af.dm, 1, function(x) sum(x > 50)))
Sys.time() - start_time


#~~~~~~~~~~~~~~~~~~~~~~#
# Generate a blocklist #
#~~~~~~~~~~~~~~~~~~~~~~#

# Some variants were also found in the K562/BT142 cell line mixing experiment and therefore possibly artifacts.
# Get mean allele frequence in K562 and BT142 cells and add it to vars.tib
backgroundvars.tib <- read_tsv("../3_Cell_line_mixes_variants/TenX_CellLineMix_Variants1.txt")
backgroundvars.tib <- backgroundvars.tib %>% select(var, mean_cov.unionCells, mean_af.K562, mean_af.BT142)

# Function to determine the VAF fold change between the two cell lines
var_fc <- function(a, b) {
    # Distance to homoplasmy
    x <- min(a, abs(100-a))
    y <- min(b, abs(100-b))
    d <- abs(a-b)
    # Fold difference
    z <- (min(x,y)+d) / min(x,y)
    return(z)
}

backgroundvars.tib$FoldChange <- apply(backgroundvars.tib, 1, function(x) var_fc(as.numeric(x[3]), as.numeric(x[4])) )

# Make a blocklist containing variants that had good coverage in K562/BT142 cells and were detected in both cell lines,
# but did not have a large difference (FoldChange) as you would expect for true variants
blocklist.var <- backgroundvars.tib %>% filter(mean_cov.unionCells > 5) %>%
    filter(between(mean_af.K562, 0.1, 99.9), between(mean_af.BT142, 0.1, 99.9), FoldChange < 5) %>% .$var
# Also block variants of position 3107, which is N in the reference
blocklist.var <- c(blocklist.var, backgroundvars.tib[grepl("N", backgroundvars.tib$var),]$var)

# Check which variants are in the blocklist visually ----------------------------------------------

# Polygon coordinates
tib <- tibble(x = c(0.1, 0.1, 80/4.8,  99.5, 99.9, 99.9, 5*100/6,  0.5, 0.1),
              y = c(0.1, 0.5, 5*100/6, 99.9, 99.9, 99.5, 80/4.8,   0.1, 0.1))
# One symbol per variant
points.tib <- backgroundvars.tib %>% filter(mean_cov.unionCells > 5) %>%
    mutate(block = var %in% blocklist.var) %>% arrange(block)

pdf("4.2_1_Blocklist.pdf")
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Test different variant selection thresholds #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Specify variant selection thresholds to test. voi = variant of interest
conditions.tib <- tibble(min_clone_size = rep(1:25, 4),
                         min_vaf = rep(c("n1", "n5", "n10", "n50"), each = 25),
                         vois = NA,
                         n_vois = NA,
                         cells = NA,
                         n_cells = NA,
                         transitions = NA)
vois.ls <- vector(mode = "list", length = nrow(conditions.tib))
cells.ls <- vector(mode = "list", length = nrow(conditions.tib))

# Quality threshold, including selection of variants that are absent from 90% of cells
vars_filter.tib <- vars.tib %>% filter(mean_cov > 5, quality >= 30, n0 > 0.9*ncol(af.dm))

# Fill in conditions.tib
for (x in 1:nrow(conditions.tib)) {
    print(x)
    
    min_clone_size <- conditions.tib$min_clone_size[x]
    min_vaf <- conditions.tib$min_vaf[x]

    # Define variants of which the number of cells exceeding min_vaf is higher than min_clone_size
    voi.ch <- vars_filter.tib$var[vars_filter.tib[[min_vaf]] >= min_clone_size]
    
    # Exclude bad variants (as well as 1583_A>G based on downstream analysis)
    voi.ch <- setdiff(voi.ch, c("1583_A>G", blocklist.var))

    # Which cells are positive for at least one of the variants?
    af_subset.dm <- af.dm[voi.ch,]
    positive_cells <- colnames( af_subset.dm[,colSums(af_subset.dm > 1)] )
    
    # Add information to summary table
    conditions.tib[x,"n_vois"] <- length(voi.ch)
    conditions.tib[x,"n_cells"] <- length(positive_cells)
    # Transitions vs. transversions
    conditions.tib[x,"transitions"] <- mean( str_count(voi.ch, "G>A|A>G|C>T|T>C") )
    
    # Save variants and cells for these parameters
    vois.ls[[x]] <- voi.ch
    cells.ls[[x]] <- positive_cells
}
conditions.tib$vois <- vois.ls
conditions.tib$cells <- cells.ls

# Visualize
pdf("4.2_2_Recovered_cells.pdf", width = 5, height = 5)
conditions.tib %>%
    mutate(min_vaf = str_c(gsub("n", "", min_vaf), "%")) %>%
    mutate(min_vaf = factor(min_vaf, levels = c("1%", "5%", "10%", "50%"))) %>%
    ggplot(aes(x = min_clone_size, y = n_cells, color = min_vaf)) + #, size = n_vois
    geom_hline(yintercept = ncol(af.dm)) +
    geom_point() +
    coord_cartesian(ylim = c(0, ncol(af.dm))) +
    ylab("Number of cells with informative variant") +
    xlab("Minimum clone size") +
    theme_bw() +
    guides(color = guide_legend(override.aes = list(size = 3),
                                title = "Minimum VAF")) +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black"))
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Deeper analysis of four variant selection thresholds #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Here we look into more detail into four parameter combinations:
# 1. Select variants present in at least 5 cells with a VAF of >10%
# 2. Select variants present in at least 10 cells with a VAF of >10%
# 3. Select variants present in at least 5 cells with a VAF of >50%
# 4. Select variants present in at least 10 cells with a VAF of >50%
conditions_subset.tib <- conditions.tib %>% filter(min_clone_size %in% c(5,10), min_vaf %in% c("n10","n50"))
conditions_subset.tib$ARI <- NA

# Pick one of the threshold settings
a <- 1
a <- 2
a <- 3
a <- 4 # this is the one we ended up using for the paper
# Or run all four
for (a in 1:nrow(conditions_subset.tib)) {
print(a)
voi.ch <- conditions_subset.tib$vois[[a]]

# List cell IDs that are positive for each of the vois --------------------------------------------
positive_cells.ls <- list()
for (v in voi.ch) {
    # Determine cells with an appreciable VAF
    current_cells.ch <- colnames(af.dm)[af.dm[v,]>1]
    # Save cell IDs for positive cells
    positive_cells.ls[[v]] <- current_cells.ch
}

# Make a tibble of cells marked by each voi
positive_cells.tib <- as_tibble(bind_rows(lapply(positive_cells.ls, function(x) data.frame(cell = x)), .id = "variant")[,2:1]) %>%
    mutate(variant = factor(variant, levels = voi.ch))
cov.mat <- as.matrix( assays(maegatk.rse)$coverage )
positive_cells.tib$cov <- apply(positive_cells.tib, 1, function(x) { cov.mat[as.numeric(cutf(x[2], d = "_")),x[1]] } )

# Add cell type and TRB columns. The TRB (T-cell receptor) sequence allows us to calculate the adjusted rand index (ARI) later
positive_cells.tib <- positive_cells.tib %>%
    left_join(dplyr::select(as_tibble(seu@meta.data, rownames = "cell"), cell, CellType, TRB_CDR3))

# Unlike the full analysis of the clonal hematopoiesis sample (4.3_Variants_Of_Interest.R), I am not grouping
# similar variants together into groups (clones), but treating each variant separately.


# Heatmap -----------------------------------------------------------------------------------------

# Prepare matrix of variants of interest in cells that are positive for at least one
af_voi.mat <- af.dm[voi.ch,]
af_subset.mat <- af_voi.mat[,apply(af_voi.mat, 2, function(x) sum(x > 1) > 0)]

# Customize column order. This is different from the strategy for K562 subclones.
plot_order.mat <- af_subset.mat
for (x in rev(voi.ch)) { plot_order.mat <- plot_order.mat[,order(-plot_order.mat[x,])] }

# Set upper VAF limit for visualization clarity
if (a %in% c(1, 2)) { plot_order.mat[plot_order.mat > 20] <- 20 }

# Add annotation bar
anno.tib <- tibble(cell = colnames(plot_order.mat)) %>% left_join(metadata.tib, by = "cell") %>% select(CellType)
ha <- HeatmapAnnotation(df = data.frame(anno.tib), col = list(CellType = mycol.ch))

# Plot
hm <- Heatmap(plot_order.mat,
              col = colorRamp2(seq(0, round(max(plot_order.mat)), length.out = 9),
                               c("#FCFCFC","#FFEDB0","#FFDF5F","#FEC510","#FA8E24","#F14C2B","#DA2828","#BE2222","#A31D1D")),
              show_row_names = ifelse(nrow(plot_order.mat) < 100, T, F),
              show_column_names = F,
              cluster_columns = F,
              cluster_rows = F,
              row_names_gp = gpar(fontsize = 10),
              name = "AF",
              heatmap_legend_param = list(border = "#000000", grid_height = unit(10, "mm")),
              top_annotation = ha,
              border = T,
              width = unit(200, "mm"),
              height = unit(100, "mm"),
              use_raster = T,
              raster_quality = 5)
pdf(str_c("4.2_3_", a ,"_Heatmap.pdf"), width = 12, height = 6)
print(hm)
dev.off()


# Calculate ARI with TRB data ---------------------------------------------------------------------

# First, apply filters similar to 5.2_Clonal_overlap.R. Some cells have multiple MT variants: only the
# variant with the highest coverage is kept to facilitate downstream processing.
filtered_cells.tib <- positive_cells.tib %>% group_by(cell) %>% 
    group_modify(~ {.x %>% arrange(desc(cov))}) %>%
    filter(!duplicated(cell))

# Determine largest mtDNA clone in TRB clones and vice versa, and the size of clones defined by both.
filtered_cells.tib <- filtered_cells.tib %>%
    ungroup %>% group_by(TRB_CDR3) %>%
    group_modify(~ {.x %>% mutate(largest_mt = max(table(variant)))}) %>%
    ungroup %>% group_by(variant) %>%
    group_modify(~ {.x %>% mutate(largest_trb = max(table(TRB_CDR3)))}) %>%
    ungroup
# Select only TRB clones that have an MT clone of >5 cells and vice versa
filtered_cells.tib <- filtered_cells.tib %>% filter(largest_mt > 5, largest_trb > 5)

conditions_subset.tib$ARI[a] <- adjustedRandIndex(filtered_cells.tib$variant, filtered_cells.tib$TRB_CDR3)

}

# Save for Supplemental Figure 11b
write_tsv(select(conditions_subset.tib, min_clone_size, min_vaf, n_vois, n_cells, transitions, ARI), file = "4.2_Four_treshold_results.txt")

# For the paper: "We plotted the largest and most distinct 23 clones using 26 informative mtDNA variants (14.9% of cells were assigned to these clones, Figure 2B)."
# Save for the next script. Make sure this section was run with a == 4 last.
write.table(voi.ch, file = "4.2_vois.txt", quote = F, col.names = F, row.names = F)


