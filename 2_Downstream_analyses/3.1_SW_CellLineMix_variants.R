# Peter van Galen, 201101
# Distinguish K562 and BT142 cells in cell line mixing experiments


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
library(readxl)
library(ggrastr)

setwd("~/DropboxPartners/Projects/Maester/AnalysisPeter/3_Cell_line_mixes_variants")
rm(list=ls())

# Functions and colors (available at https://github.com/vangalenlab/MAESTER-2021)
source("../210215_FunctionsGeneral.R")
popcol.df <- read_excel("../MAESTER_colors.xlsx")
heatcol.ch <- read_excel("../MAESTER_colors.xlsx", sheet = 2, col_names = "heatcol")$heatcol

# Import data for Seq-Well S^3 Cell Line Mix (available at https://vangalenlab.bwh.harvard.edu/maester-2021/)
seu <- readRDS("../2_Cell_line_mixes_decontX/SW_CellLineMix_Seurat_Keep.rds")
se.ls <- readRDS("../1_MT_Coverage/SW_CellLineMix_mr3_maegatk.rds")


#~~~~~~~~~~~~~~~~~~~~#
#### Prepare data ####
#~~~~~~~~~~~~~~~~~~~~#

# Combine information
names(se.ls)
common.cells <- intersect(colnames(seu), colnames(se.ls[[2]]))

# Fraction of scRNA-seq cells captured with Maester (99.4%):
rnaseq_cell_number <- ncol(seu)
common_cell_number <- length(common.cells)
common_cell_number / rnaseq_cell_number

# Put all objects in the same order
seu <- seu[,common.cells]
maegatk.rse <- se.ls[[2]][,common.cells]

# Prepare allele frequency matrix. Rows represents a position along the mitochondrial genome and the three possible disagreements with the reference (except 3107 has four possible disagreements because the reference is N)
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))*100

# Collect cell information 
cells.tib <- tibble(cell = common.cells,
                    orig.ident = seu$orig.ident,
                    CellType_RNA = seu$CellType,
                    UMAP_1 = seu$UMAP_1,
                    UMAP_2 = seu$UMAP_2,
                    Mean_Cov = maegatk.rse$depth)

# Make groups of cell IDs
CellSubsets.ls <- list(unionCells = cells.tib$cell,
                       K562 = filter(cells.tib, CellType_RNA == "K562")$cell,
                       BT142 = filter(cells.tib, CellType_RNA == "BT142")$cell)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Collect variant information ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# This takes a while. Consider skipping to read_tsv at the end.

# Get the mean allele frequency and coverage for every cell subset.
mean_af.ls <- lapply(CellSubsets.ls, function(x) rowMeans(af.dm[,x]))
mean_cov.ls <- lapply(CellSubsets.ls, function(x) rowMeans(assays(se.ls[[2]])[["coverage"]][,x])[as.numeric(cutf(rownames(af.dm), d = "_"))])
names(mean_af.ls) <- paste0("mean_af.", names(mean_af.ls))
names(mean_cov.ls) <- paste0("mean_cov.", names(mean_cov.ls))

# Get the quantiles of the VAFs of each variant in each cell subset
quantiles <- c("q01" = 0.01, "q10" = 0.1, "q50" = 0.5, "q90" = 0.9, "q99" = 0.99)
start_time <- Sys.time()
quantiles.ls <- lapply(quantiles, function(x) lapply(CellSubsets.ls, function(y) apply(af.dm[,y], 1, quantile, x) ))
Sys.time() - start_time

# Get the mean quality for each variant. This takes an hour.
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
write_tsv(vars.tib, "SW_CellLineMix_Variants1.txt")
vars.tib <- read_tsv("SW_CellLineMix_Variants1.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Detect homoplasmic variants ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Variants 20-fold coverage, VAF of <10% in the bottom 10% of cells and VAF of >90% in top 10% of cells
voi.ch <- vars.tib %>% filter(mean_cov.unionCells > 20,
                              quality > 30,
                              q10.unionCells < 10,
                              q90.unionCells > 90) %>% .$var
# Assess transitions vs. transversions
str_view(voi.ch, "G>A|A>G|C>T|T>C"); mean( str_count(voi.ch, "G>A|A>G|C>T|T>C") )

# Generate matrices with coverage, allele frequency and reference / variant reads
cov_voi.mat <- assays(maegatk.rse)[["coverage"]][as.numeric(cutf(voi.ch, d = "_")),]
af_voi.mat <- af.dm[voi.ch,]

### Plot VAF in each cell, sorted from low to high, illustrating the selection process
pdf(paste0(seu@project.name, "_1_Sorted_VAFs.pdf"), width = 10, height = 5)
par(mar=c(4,4,2,8), xpd = T)
plot(NA, xlim = c(0, ncol(af_voi.mat)), ylim = c(0, 100), xlab = "Cells (sorted separately for each variant)", ylab = "Variant Allele Frequency")
for (n in 1:length(voi.ch)) {
    v <- voi.ch[n]
    points(as.numeric(sort(af_voi.mat[v,])), pch = 16, col = popcol.df$hex[n])
    text(x = 1.1*ncol(af_voi.mat), y = 110-10*n, label = v, col = popcol.df$hex[n], pos = 4)
}
dev.off()

# Add coverage and allele frequency info from variants of interest to cells.tib
for (voi in voi.ch) {
cells.tib <- cells.tib %>%
    left_join(as_tibble(assays(maegatk.rse)[["coverage"]][as.numeric(cutf(voi, d = "_")),], rownames = "cell"), by = "cell") %>%
    left_join(as_tibble(af.dm[voi,], rownames = "cell"), by = "cell") %>%
    rename(value.x = str_c("cov_", str_replace(voi, ">", ".")), value.y = str_c("af_", str_replace(voi, ">", ".")))
}

### For each variant of interest, plot UMAP of cells, colored by VAF
pdf(paste0(seu@project.name, "_2_VAF_UMAP.pdf"))
for (voi in voi.ch) {
    #voi <- voi.ch[1]
    cov_colname <- str_c("cov_", str_replace(voi, ">", "."))
    af_colname <- str_c("af_", str_replace(voi, ">", "."))
    
    # Select cells with three genotyped transcripts, then plot
    print(
    cells.tib %>% filter(.[[cov_colname]] > 3) %>%
        ggplot(aes_string(x = "UMAP_1", y = "UMAP_2", color = af_colname)) +
        geom_point_rast() +
        scale_color_gradientn(colors = heatcol.ch[2:10]) +
        theme_classic() + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
        ggtitle(voi)
    )
}
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Classify cells by homoplasmic variants ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Function to sum up the reads supporting either reference or mutant allele for a particulal variant of interest
# voi is the variant of interest for which to sum up all the supporting reads (e.g. "709_G>A" or "709_G.A")
# maegatk is the output from the maegatk software of class RangedSummarizedExperiment
supporting_calls <- function(voi, maegatk) {
    #voi <- voi.ch[1]
    #maegatk <- maegatk.rse
    # From the variant, determine position, reference and mutant allele
    pos <- as.numeric(cutf(voi, d = "_"))
    ref <- cutf(voi, d = "_|>|\\.", f = 2)
    mut <- cutf(voi, d = ">|\\.", f = 2)
    
    # Put all the counts for this variant from the maegatk object in a data frame
    counts.df <- do.call(cbind, lapply(assays(maegatk)[grepl("counts", names(assays(maegatk)))], function(x) x[pos,]))
    
    # Add up reads supporting reference or mutated allele
    support.df <- data.frame(ref = rowSums(counts.df[,grepl(paste0("^", ref), colnames(counts.df))]),
                             mut = rowSums(counts.df[,grepl(paste0("^", mut), colnames(counts.df))]))
    
    support.df
}

# Calculate number of supporting reads for each variant; then combine
supporting_calls.ls <- lapply(voi.ch, function(x) supporting_calls(x, maegatk.rse))
names(supporting_calls.ls) <- voi.ch
combined_supporting_calls.df <- do.call(cbind, supporting_calls.ls)

# Add up supporting reads: is it a K562 or BT142 cell?
supporting_calls.tib <- tibble(cell = rownames(combined_supporting_calls.df),
    K562_supporting_calls = rowSums(combined_supporting_calls.df[,c("709_G>A.mut", "1888_G>A.mut", "1420_T>C.mut", "2141_T>C.mut", "9117_T>C.mut", "7990_C>T.ref")]),
    BT142_supporting_calls = rowSums(combined_supporting_calls.df[,c("709_G>A.ref", "1888_G>A.ref", "1420_T>C.ref", "2141_T>C.ref", "9117_T>C.ref", "7990_C>T.mut")]))
           
# Combine information
cells.tib <- left_join(cells.tib, supporting_calls.tib, by = "cell")

# Label cells as donor, host, or contaminated based on the number of supporting reads
#cells.tib %>% mutate(reads = K562_supporting_calls + BT142_supporting_calls) %>% arrange(reads)
cells.tib <- cells.tib %>% mutate(CellType_MT = ifelse(K562_supporting_calls + BT142_supporting_calls < 30, yes = "NoCoverage", no =
    ifelse(K562_supporting_calls > 3 & K562_supporting_calls > 10*BT142_supporting_calls, yes = "K562", no =
    ifelse(BT142_supporting_calls > 3 & BT142_supporting_calls > 10*K562_supporting_calls, yes = "BT142", no = "Contaminated")))) %>%
    mutate(CellType_MT = factor(CellType_MT, levels = c("Contaminated", "BT142", "K562", "NoCoverage")))

# Plot supporting reads
pdf(paste0(seu@project.name, "_3_SupportingReads.pdf"))
cells.tib %>% arrange(CellType_MT) %>%
    ggplot(aes(x = K562_supporting_calls, y = BT142_supporting_calls, color = CellType_MT)) +
    coord_cartesian(xlim = c(0,4000)) +
    geom_point_rast(alpha = 1, size = 2) +
    scale_color_manual(values = c("#E17973", "#8E87F5", "#7BF581", NA)) +
    theme_classic() +
    theme(aspect.ratio = 1)
dev.off()

### Pie chart of classifications
table(cells.tib$CellType_MT, useNA = "ifany")
frac.tib <- tribble(~CellType_MT, ~color, ~value,
                    "Contaminated", "#E17973", sum(cells.tib$CellType_MT == "Contaminated"),
                    "BT142", "#8E87F5", sum(cells.tib$CellType_MT == "BT142"),
                    "K562", "#7BF581", sum(cells.tib$CellType_MT == "K562"),
                    "NoCoverage", "black", sum(cells.tib$CellType_MT == "NoCoverage") + rnaseq_cell_number - common_cell_number)
frac.tib <- frac.tib %>% mutate(CellType_MT = factor(CellType_MT, levels = c("Contaminated", "BT142", "K562", "NoCoverage")))
sum(frac.tib$value) == rnaseq_cell_number # Ok

# Plot
pdf(paste0(seu@project.name, "_4_Piechart.pdf"))
ggplot(frac.tib, aes(x = "", y = value, fill = CellType_MT)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    scale_fill_manual(values = frac.tib$color, labels = paste0(frac.tib$CellType_MT, " (", frac.tib$value, ")")) +
    theme_void()
dev.off()


### UMAP with CellType_RNA from RNA-seq and CellType_MT from MT DNA
pdf(paste0(seu@project.name, "_5_UMAP_classify.pdf"))
print(
cells.tib %>% filter(! CellType_MT %in% c("Contaminated", "NoCoverage")) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = CellType_MT)) +
    geom_point_rast() +
    scale_color_manual(values = c("#575757", "#ff9233")) +
    theme_classic() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    ggtitle("Classification from mitochondrial variants")
)
print(
cells.tib %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = CellType_RNA)) +
    geom_point_rast() +
    scale_color_manual(values = c("#575757", "#ff9233")) +
    theme_classic() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    ggtitle("Cell type from RNA-seq")
)
dev.off()


#~~~~~~~~~~~~~~~~#
#### Heatmaps ####
#~~~~~~~~~~~~~~~~#

# What cells were called and have coverage of three transcripts at each of the variants?
cells.ch <- cells.tib %>% filter(CellType_MT %in% c("K562", "BT142")) %>%
    select(cell, starts_with("cov_")) %>%
    filter(apply(.[,-1], 1, function(x) all(x > 3))) %>%
    .$cell

# Create heatmap annotation (RNA-seq clusters)
ha <- HeatmapAnnotation(CellType_RNA = cells.tib$CellType_RNA[cells.ch],
                        col = list(CellType_RNA = c("K562" = "#7BF581", "BT142" = "#8E87F5")))

# Visualize: clustered VAFs
hm1 <- Heatmap(af_voi.mat[,cells.ch],
               col=c("#FCFCFC","#FFEDB0","#FFDF5F","#FEC510","#FA8E24","#F14C2B","#DA2828","#BE2222","#A31D1D"),
               show_row_names = T,
               show_column_names = F,
               cluster_columns = T,
               cluster_rows = F,
               name = "AF",
               top_annotation = ha,
               border = T,
               width = unit(3.5, "inch"),
               height = unit(3.5, "inch"),
               heatmap_legend_param = list(border = "black", at = c(0,25,50,75,100)),
               use_raster = T)

# Correlate cells based on informative variants
cor.mat <- cor(af_voi.mat[,cells.ch])
cor.mat[is.na(cor.mat)] <- 1

hm2 <- Heatmap(cor.mat,
               col=c("#FCFCFC","#FFEDB0","#FFDF5F","#FEC510","#FA8E24","#F14C2B","#DA2828","#BE2222","#A31D1D"),
               show_row_names = F, 
               show_column_names = F,
               name = "MitoCor",
               top_annotation = ha,
               border = T,
               width = unit(3.5, "inch"),
               height = unit(3.5, "inch"),
               heatmap_legend_param = list(border = "black", at = c(-1,0,1)),
               use_raster = T)

pdf(str_c(seu@project.name, "_6_Heatmaps.pdf"))
print(hm1)
print(hm2)
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~#
#### Metrics & save ####
#~~~~~~~~~~~~~~~~~~~~~~#

# Agreement between CellType_RNA and CellType_MT
cells.tib$CellType_RNA %>% table(useNA= "always")
cells.tib$CellType_MT %>% table(useNA = "always")
cells.tib %>% filter(CellType_MT == "K562") %>% .$CellType_RNA %>% table(useNA = "always")
cells.tib %>% filter(CellType_MT == "BT142") %>% .$CellType_RNA %>% table(useNA = "always")
# "Of the cells that were classified by mitochondrial variants, 2,129/2,130 were concordant with mRNA clusters for Seq-Well"

# Remove unneccessary information
cells.tib <- cells.tib[,!str_detect(colnames(cells.tib), str_replace_all(str_c(c(voi.ch, "supporting"), collapse = "|"), ">", "."))]
write_rds(cells.tib, "SW_CellLineMix_Cells.rds")


