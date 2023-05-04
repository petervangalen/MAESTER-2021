# Peter van Galen, 210704
# Compare MAESTER and bulk ATAC-seq variants in K562 cells

# Prerequisites
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(tidyverse)
library(ggrepel)
library(gridExtra)

setwd("~/DropboxMGB/Projects/Maester/AnalysisPeter/9_K562_ATAC/")

rm(list=ls())

# Functions (available at https://github.com/vangalenlab/MAESTER-2021)
source("../Auxiliary_files/210215_FunctionsGeneral.R")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Import and wrangle bam-readount results #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# bam-readcount was used to call variants in the K562 bulk ATAC-seq data. Load the results. 
counts.df <- read.table("bam-readcount_results.txt", sep = "\t", fill = T, col.names = c("chr", "pos", "ref", "depth", "base", "A", "C", "G", "T", "N", "Other1", "Other2"))
# Importing these results is slightly unusual so you should check that it all went well
nrow(counts.df) == 16569
cutf(counts.df$A, d = ":") %>% table
cutf(counts.df$C, d = ":") %>% table
cutf(counts.df$G, d = ":") %>% table
cutf(counts.df$`T`, d = ":") %>% table
cutf(counts.df$N, d = ":") %>% table
cutf(counts.df$Other1, d = ":") %>% table
cutf(counts.df$Other2, d = ":") %>% table


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Get variants from SW and 10X MAESTER data  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

### Seq-Well (the scripts to generate these files are available on Github)
all_sw.tib <- read_tsv("../3_Cell_line_mixes_variants/SW_CellLineMix_Variants1.txt")
all_sw.tib <- all_sw.tib %>% select(var, mean_cov.unionCells, mean_af.K562) %>% dplyr::rename(SW_coverage = mean_cov.unionCells, SW_VAF = mean_af.K562)
# Get variants shown in Figure 1F and H of the first submission
homoplasmic_sw.tib <- read_tsv("../3_Cell_line_mixes_variants/Figure1F_vars.txt")
subclonal_sw.tib <- read_tsv("../3_Cell_line_mixes_variants/Figure1H_vars.txt")

### 10x
all_10x.tib <- read_tsv("../3_Cell_line_mixes_variants/TenX_CellLineMix_Variants1.txt")
all_10x.tib <- all_10x.tib %>% select(var, mean_cov.unionCells, mean_af.K562) %>% dplyr::rename(TenX_coverage = mean_cov.unionCells, TenX_VAF = mean_af.K562)
# Get variants shown in Figure S7 of the first submission
homoplasmic_10x.tib <- read_tsv("../3_Cell_line_mixes_variants/FigureS7C_vars.txt")
subclonal_10x.tib <- read_tsv("../3_Cell_line_mixes_variants/FigureS7E_vars.txt")

# Combine
maester.tib <- tibble(var = sort( unique( c(homoplasmic_sw.tib$var, subclonal_sw.tib$var,
                                            homoplasmic_10x.tib$var, subclonal_10x.tib$var) ) ) )
maester.tib <- maester.tib %>% left_join(all_sw.tib, by = "var") %>% left_join(all_10x.tib, by = "var")
# Add some informative columns
maester.tib <- maester.tib %>%
    mutate(pos = as.numeric(cutf(var, d = "_", f = 1)), ref = cutf(var, d = "_|>", f = 2), alt = cutf(var, d = ">", f = 2)) %>%
    filter(var != "7990_C>T") # this one is homoplasmic in BT142


#~~~~~~~~~~~~~~~~~~#
# Combine and plot #
#~~~~~~~~~~~~~~~~~~#

# Fill a matrix of variants from MAESTER with counts from ATAC (bam-readcount)
counts.mat <- matrix(nrow = nrow(maester.tib), ncol = 4)
colnames(counts.mat) <- c("pos", "depth", "ref_count", "alt_count")
for (n in 1:nrow(maester.tib)) {
    #n <- 1
    pos <- maester.tib[n,]$pos
    ref <- maester.tib[n,]$ref
    alt <- maester.tib[n,]$alt
    
    stopifnot(counts.df[pos,]$pos == pos, counts.df[pos,]$ref == ref)
    
    counts.mat[n, "pos"] <- pos
    counts.mat[n, "depth"] <- counts.df[pos,]$depth
    counts.mat[n, "ref_count"] <- as.numeric( cutf(counts.df[pos,ref], d = ":", f = 2) )
    counts.mat[n, "alt_count"] <- as.numeric( cutf(counts.df[pos,alt], d = ":", f = 2) )
}

# Merge information
merge.tib <- maester.tib %>% left_join(as_tibble(counts.mat), by = "pos")
merge.tib <- merge.tib %>% mutate(ATAC_VAF = alt_count/depth*100)

# Save data
final_summary.tib <- merge.tib %>% select(var, ATAC_VAF, SW_VAF, TenX_VAF, depth, SW_coverage, TenX_coverage) %>%
    dplyr::rename(ATAC_coverage = depth) %>%
    arrange(desc(ATAC_VAF)) %>%
    mutate(index = row_number(), .before = 1)
write_tsv(final_summary.tib, file = "final_summary.txt")

# Based on this comparison, 5378A>G, 6384G>A, 7693C>T and 8251G>A (from the bioRxiv preprint) are removed in the scripts 3.2_SW_K562_clones.R and 3.4_TenX_K562_clones.R

# Plot
g1 <- merge.tib %>%
    ggplot(aes(x = SW_VAF, y = ATAC_VAF, label = gsub("_", "", var))) +
    #geom_text_repel(color = "#66cdaa", size = 4, segment.size = 0.2) +
    geom_point() +
    scale_x_continuous(trans = "log10", limits = c(0.1, 100)) +
    scale_y_continuous(trans = "log10", limits = c(0.1, 100)) +
    xlab("Seq-Well MAESTER VAF") +
    ylab("Bulk ATAC VAF") +
    ggtitle(str_c("R = ", round(cor(merge.tib$SW_VAF, merge.tib$ATAC_VAF), 4))) +
    theme_bw() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))

g2 <- merge.tib %>%
    ggplot(aes(x = TenX_VAF, y = ATAC_VAF, label = gsub("_", "", var))) +
    #geom_text_repel(color = "#66cdaa", size = 4, segment.size = 0.2) +
    geom_point() +
    scale_x_continuous(trans = "log10", limits = c(0.1, 100)) +
    scale_y_continuous(trans = "log10", limits = c(0.1, 100)) +
    xlab("10x MAESTER VAF") +
    ylab("Bulk ATAC VAF") +
    ggtitle(str_c("R = ", round(cor(merge.tib$ATAC_VAF, merge.tib$TenX_VAF), 4))) +
    theme_bw() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))

g3 <- merge.tib %>%
    ggplot(aes(x = SW_VAF, y = TenX_VAF, label = gsub("_", "", var))) +
    #geom_text_repel(color = "#66cdaa", size = 4, segment.size = 0.2) +
    geom_point() +
    scale_x_continuous(trans = "log10", limits = c(0.1, 100)) +
    scale_y_continuous(trans = "log10", limits = c(0.1, 100)) +
    xlab("Seq-Well MAESTER VAF") +
    ylab("10x MAESTER VAF") +
    ggtitle(str_c("R = ", round(cor(merge.tib$SW_VAF, merge.tib$TenX_VAF), 4))) +
    theme_bw() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))

pdf("9.3_Correlations.pdf", width = 10, height = 4)
grid.arrange(g1, g2, g3, nrow = 1)
dev.off()



