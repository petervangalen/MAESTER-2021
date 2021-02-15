# Daniel Ssozi, 201214
# Create Seurat objects for all samples with all genotyping information
# In a way that is easy to expand if we get new data
# For BPDCN628, BPDCN712 and BPDCN712R

options(stringsAsFactors = F)
options(scipen = 999)

library(tidyverse)
library(gdata)
library(Seurat)
library(KernSmooth)
library(data.table)

rm(list=ls())

# Set working directory
setwd("~/DropboxPartners/Projects/Single-cell_BPDCN/AnalysisDaniel/") ### For Peter
setwd("/Users/dz855/Dropbox (Partners HealthCare)/Single-cell_BPDCN/AnalysisDaniel") ### For Daniel
dir.ch <- "201214_Seurat_GoT/"
dir.create(dir.ch)

# Load commonly used function
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Load Seurat object
BPDCN628 <- readRDS("../AnalysisPeter/200508_RandomForest/BPDCN628_Seurat_Predict.rds")
BPDCN712 <- readRDS("../AnalysisPeter/200508_RandomForest/BPDCN712_Seurat_Predict.rds")
NormalRelapse_old <- readRDS("../AnalysisPeter/200726_RelapseViz/NormalRelapse_Seurat.rds")
# Need to split for NormalRelapse cells so that duplicate cell barcodes between BPDCN712R and BM are not automatically overwritten with "no call" in the annotate_Seurat_with_GoT function
BM <- subset(NormalRelapse_old, orig.ident == "BM")
BPDCN712R <- subset(NormalRelapse_old, orig.ident == "BPDCN712R")


################################
# Load results from GoT pipeline
################################

# Where to find data + amino acid change
FilteredCells_files.tib <- tribble(~Sample, ~FilteredCells, ~Mut, ~Note,
    ### BPDCN628
    "BPDCN628", "PCRanalysis/200620_BPDCN628/CUX1.2707/CUX1.2707.FilteredCells.txt", "CUX1.L911fs*", "CH mutation from RHP", 
    "BPDCN628", "PCRanalysis/200620_BPDCN628/EZH2.1490/EZH2.1490.FilteredCells.txt", "EZH2.A505T", "Germline mutation",
    "BPDCN628", "PCRanalysis/200620_BPDCN628/TET2.4291/TET2.4291.FilteredCells.txt", "TET2.E1437fs*", "CH mutation from RHP", 
    "BPDCN628", "PCRanalysis/200620_BPDCN628/TET2.4592/TET2.4592.FilteredCells.txt", "TET2.Q1547*", "CH mutation from RHP", 
    
    ### BPDCN712 dx from Benchling note: "BPDCN712 genotyping (RHP mutations)"
    #"BPDCN712", "PCRanalysis/200108_BPDCN712/ASXL1.1898/ASXL1.1898.FilteredCells.txt", "ASXL1.G642fs.1", # Do not use for paper (use primer ASXL1-1886 instead; similar results)
    "BPDCN712", "PCRanalysis/200108_BPDCN712/TET2.2340/TET2.2340.FilteredCells.txt", "TET2.S792*", "CH mutation from RHP", 
    "BPDCN712", "PCRanalysis/200108_BPDCN712/TET2.3078/TET2.3078.FilteredCells.txt", "TET2.Q1034*", "CH mutation from RHP",
    "BPDCN712", "PCRanalysis/200108_BPDCN712/TET2.3626/TET2.3626.FilteredCells.txt", "TET2.R1216*", "CH mutation from RHP",
    "BPDCN712", "PCRanalysis/200108_BPDCN712/TET2.4104/TET2.4104.FilteredCells.txt", "TET2.H1380Y", "CH mutation from RHP",
    
    ### BPDCN712 dx from "BPDCN712 and BPDCN712R genotyping (comprehensive)"
    "BPDCN712", "PCRanalysis/201127_BPDCN712/ACAP2.239/ACAP2.239.FilteredCells.txt", "ACAP2.L97M", "Skin-specific mutation",
    "BPDCN712", "PCRanalysis/201127_BPDCN712/ASXL1.1886/ASXL1.1886.FilteredCells.txt", "ASXL1.G642fs", "CH mutation from RHP",
    "BPDCN712", "PCRanalysis/201127_BPDCN712/CWF19L2.613/CWF19L2.613.FilteredCells.txt", "CWF19L2.K221R", "Skin-specific mutation", # Exclude from some analyses b/c not UV
    "BPDCN712", "PCRanalysis/201127_BPDCN712/DOLPP1.633/DOLPP1.633.FilteredCells.txt", "DOLPP1.R227S", "Skin-specific mutation", # Exclude from some analyses b/c not UV
    "BPDCN712", "PCRanalysis/201127_BPDCN712/HNRNPUL1.1629/HNRNPUL1.1629.FilteredCells.txt", "HNRNPUL1.F559F", "Skin-specific mutation",
    #"PCRanalysis/201127_BPDCN712/IKZF1.1477/IKZF1.1477.FilteredCells.txt", "IKZF1.S507L", # IronThrone-GoT SummTable not generated for unclear reasons
    "BPDCN712", "PCRanalysis/201127_BPDCN712/MAP4K5.1949/MAP4K5.1949.FilteredCells.txt", "MAP4K5.P667S", "Skin-specific mutation",
    "BPDCN712", "PCRanalysis/201219_MTAP/MTAP.308.1/MTAP.308.1.FilteredCells.txt", "MTAP.rearr.1", "Skin-specific mutation indicating CDKN2A deletion",
    "BPDCN712", "PCRanalysis/201219_MTAP/MTAP.308.2/MTAP.308.2.FilteredCells.txt", "MTAP.rearr.2", "Skin-specific mutation indicating CDKN2A deletion",
    "BPDCN712", "PCRanalysis/201219_MTAP/MTAP.308.3/MTAP.308.3.FilteredCells.txt", "MTAP.rearr.3", "Skin-specific mutation indicating CDKN2A deletion",
    "BPDCN712", "PCRanalysis/201219_MTAP/MTAP.308.4/MTAP.308.4.FilteredCells.txt", "MTAP.rearr.4", "Skin-specific mutation indicating CDKN2A deletion",
    "BPDCN712", "PCRanalysis/201219_MTAP/MTAP.308.5/MTAP.308.5.FilteredCells.txt", "MTAP.rearr.5", "Skin-specific mutation indicating CDKN2A deletion",
    #"BPDCN712", "PCRanalysis/201127_BPDCN712/NIPA1.1/NIPA1.1.FilteredCells.txt", "NIPA1.A14G", "Skin-specific mutation", # Likely sequencing error (GC-rich region), only wt transcripts
    "BPDCN712", "PCRanalysis/201127_BPDCN712/RAB9A_1M/RAB9A.575.FilteredCells.txt", "RAB9A.3pUTR", "Skin-specific mutation on X chromosome", # Only used first 250K reads
    #"BPDCN712", "PCRanalysis/201127_BPDCN712/RBM12B.533/RBM12B.533.FilteredCells.txt", "RBM12B.G195S", "Skin-specific mutation", # Exclude from paper because relapse-specific
    #"BPDCN712", "PCRanalysis/201127_BPDCN712/ZNF750.1562/ZNF750.1562.FilteredCells.txt", "ZNF750.E529L", "Skin-specific mutation", # Didn't work (only 1 cell)
    
    ### BPDCN712R from "BPDCN712 and BPDCN712R genotyping (comprehensive)"
    "BPDCN712R", "PCRanalysis/201124_BPDCN712R/ACAP2.239/ACAP2.239.R.FilteredCells.txt", "ACAP2.L97M", "Skin-specific mutation",
    "BPDCN712R", "PCRanalysis/201124_BPDCN712R/ASXL1.1886/ASXL1.1886.R.FilteredCells.txt", "ASXL1.G642fs", "CH mutation from RHP",
    "BPDCN712R", "PCRanalysis/201124_BPDCN712R/CWF19L2.613/CWF19L2.613.R.FilteredCells.txt", "CWF19L2.K221R", "Skin-specific mutation", # Exclude from some analyses b/c not UV
    "BPDCN712R", "PCRanalysis/201124_BPDCN712R/DOLPP1.633/DOLPP1.633.R.FilteredCells.txt", "DOLPP1.R227S", "Skin-specific mutation", # Exclude from some analyses b/c not UV
    "BPDCN712R", "PCRanalysis/201124_BPDCN712R/HNRNPUL1.1629/HNRNPUL1.1629.R.FilteredCells.txt", "HNRNPUL1.F559F", "Skin-specific mutation",
    #"PCRanalysis/201124_BPDCN712R/IKZF1.1477/IKZF1.1477.R.FilteredCells.txt", "IKZF1.S507L", # IronThrone-GoT SummTable not generated for unclear reasons
    "BPDCN712R", "PCRanalysis/201124_BPDCN712R/MAP4K5.1949/MAP4K5.1949.R.FilteredCells.txt", "MAP4K5.P667S", "Skin-specific mutation",
    "BPDCN712R", "PCRanalysis/201219_MTAP/MTAP.308_1M.1R/MTAP.308.1R.FilteredCells.txt", "MTAP.rearr.1", "Skin-specific mutation indicating CDKN2A deletion", # Only used first 250K reads
    "BPDCN712R", "PCRanalysis/201219_MTAP/MTAP.308_1M.2R/MTAP.308.2R.FilteredCells.txt", "MTAP.rearr.2", "Skin-specific mutation indicating CDKN2A deletion", # Only used first 250K reads
    "BPDCN712R", "PCRanalysis/201219_MTAP/MTAP.308_1M.3R/MTAP.308.3R.FilteredCells.txt", "MTAP.rearr.3", "Skin-specific mutation indicating CDKN2A deletion", # Only used first 250K reads
    "BPDCN712R", "PCRanalysis/201219_MTAP/MTAP.308_1M.4R/MTAP.308.4R.FilteredCells.txt", "MTAP.rearr.4", "Skin-specific mutation indicating CDKN2A deletion", # Only used first 250K reads
    "BPDCN712R", "PCRanalysis/201219_MTAP/MTAP.308_1M.5R/MTAP.308.5R.FilteredCells.txt", "MTAP.rearr.5", "Skin-specific mutation indicating CDKN2A deletion", # Only used first 250K reads
    #"BPDCN712R", "PCRanalysis/201124_BPDCN712R/NIPA1.1/NIPA1.1.R.FilteredCells.txt", "NIPA1.A14G", "Skin-specific mutation", # Likely sequencing error (GC-rich region), only wt transcripts
    "BPDCN712R", "PCRanalysis/201124_BPDCN712R/RAB9A_1M.R/RAB9A.575.FilteredCells.txt", "RAB9A.3pUTR", "Skin-specific mutation on X chromosome", # Only used first 250K reads
    #"BPDCN712R", "PCRanalysis/201124_BPDCN712R/RBM12B.533/RBM12B.533.R.FilteredCells.txt", "RBM12B.G195S", "Skin-specific mutation", # Exclude from paper because relapse-specific
    "BPDCN712R", "PCRanalysis/201124_BPDCN712R/TET2.2340/TET2.2340.R.FilteredCells.txt", "TET2.S792*", "CH mutation from RHP (clone 1)",
    "BPDCN712R", "PCRanalysis/201124_BPDCN712R/TET2.3078/TET2.3078.R.FilteredCells.txt", "TET2.Q1034*", "CH mutation from RHP (clone 1)",
    "BPDCN712R", "PCRanalysis/201124_BPDCN712R/TET2.3626/TET2.3626.R.FilteredCells.txt", "TET2.R1216*", "CH mutation from RHP (clone 2)",
    "BPDCN712R", "PCRanalysis/201124_BPDCN712R/TET2.4104/TET2.4104.R.FilteredCells.txt", "TET2.H1380Y", "CH mutation from RHP (clone 2)"
    #"BPDCN712R", "PCRanalysis/201124_BPDCN712R/ZNF750.1562/ZNF750.1562.R.FilteredCells.txt", "ZNF750.E529L", "Skin-specific mutation" # Didn't work (only 2 cells)
)

# Save
all( FilteredCells_files.tib$FilteredCells %>% file.exists() )
write_tsv(FilteredCells_files.tib, file = str_c(dir.ch, "FilteredCells_files.txt"))


#=========================================================================#
# Function to add genotyping data from FilteredCells.txt to Seurat object #
#=========================================================================#

annotate_Seurat_with_GoT <- function(Seurat_object, FilteredCells, Mut) {
    #Seurat_object <- BPDCN712R
    #FilteredCells <- "PCRanalysis/201124_BPDCN712R/RAB9A_1M.R/RAB9A.575.FilteredCells.txt"
    #Mut <- "RAB9A.3pUTR"
    
    # Read and summarize UMI detection
    UMIs.df <- read.table(FilteredCells, header = T)
    calls.tib <- tibble(CB = UMIs.df$BC, call = factor( ifelse(UMIs.df$mutUMIs > 0, yes = "mutant", no = "wildtype"), levels = c("mutant", "wildtype", "no call")))
    
    # Generate a tibble with the same CBs as the Seurat object and genotyping calls
    genotyping.tib <- tibble(CB = cutf(colnames(Seurat_object), d = "-", f = 1))
    genotyping.tib <- genotyping.tib %>% left_join(calls.tib)
    genotyping.tib <- genotyping.tib %>% mutate(call = replace_na(call, "no call"))
    
    # Any cell barcode that occurred in multiple libraries will show up here as a duplicated cell barcode and cannot be called because we can't know which well it came from
    dups.ch <- genotyping.tib$CB[duplicated(genotyping.tib$CB)]
    genotyping.tib$call[genotyping.tib$CB %in% dups.ch] <- "no call" 
    
    # Check that all the cells are still in the same order
    stopifnot( all( genotyping.tib$CB == cutf(colnames(Seurat_object), d = "-") ) )
    
    # Add metadata to Seurat object
    Seurat_object@meta.data[,Mut] <- genotyping.tib$call
    
    return(Seurat_object)
}

# Example execution
#BPDCN712 <- annotate_Seurat_with_GoT(BPDCN712, "PCRanalysis/200108_BPDCN712/TET2.2340/TET2.2340.FilteredCells.txt", "TET2.S792*")


#====================================================#
#### Add genotyping data to Seurat objects & save ####
#====================================================#

### Add meta data for every mutation to BPDCN628
for (x in 1:nrow(filter(FilteredCells_files.tib, Sample == "BPDCN628")) ) {
    current.tib <- filter(FilteredCells_files.tib, Sample == "BPDCN628")[x,]
    BPDCN628 <- annotate_Seurat_with_GoT(Seurat_object = get(current.tib[[1]]),
                                         FilteredCells = current.tib[[2]],
                                         Mut = current.tib[[3]])
}
BPDCN628@meta.data %>% head
saveRDS(BPDCN628, file = paste0(dir.ch, "BPDCN628_Seurat_Genotyped.rds"))


### Add meta data for every mutation to BPDCN712
for (x in 1:nrow(filter(FilteredCells_files.tib, Sample == "BPDCN712")) ) {
    current.tib <- filter(FilteredCells_files.tib, Sample == "BPDCN712")[x,]
    BPDCN712 <- annotate_Seurat_with_GoT(Seurat_object = get(current.tib[[1]]),
                                         FilteredCells = current.tib[[2]],
                                         Mut = current.tib[[3]])
}

# Special processing for MTAP:
any.mutant <- BPDCN712@meta.data %>% select(contains("MTAP")) %>% apply(1, function(x) any(grepl("mutant", x)))
# Wild type is by definition the same for all of them, and here I'm combining all mutant calls:
BPDCN712$MTAP.rearr <- BPDCN712$MTAP.rearr.1
BPDCN712$MTAP.rearr[any.mutant] <- "mutant"
# Remove separate MTAP columns
BPDCN712@meta.data[,which(grepl("MTAP.rearr.\\d", colnames(BPDCN712@meta.data)))] <- NULL

saveRDS(BPDCN712, file = paste0(dir.ch, "BPDCN712_Seurat_Genotyped.rds"))


### Add meta data for every mutation to BPDCN712R
for (x in 1:nrow(filter(FilteredCells_files.tib, Sample == "BPDCN712R")) ) {
    current.tib <- filter(FilteredCells_files.tib, Sample == "BPDCN712R")[x,]
    BPDCN712R <- annotate_Seurat_with_GoT(Seurat_object = get(current.tib[[1]]),
                                          FilteredCells = current.tib[[2]],
                                          Mut = current.tib[[3]])
}

# Special processing for MTAP
any.mutant <- BPDCN712R@meta.data %>% select(contains("MTAP")) %>% apply(1, function(x) any(grepl("mutant", x)))
# Wild type is by definition the same for all of them, and here I'm combining all mutant calls:
BPDCN712R$MTAP.rearr <- BPDCN712R$MTAP.rearr.1
BPDCN712R$MTAP.rearr[any.mutant] <- "mutant"
# Remove separate MTAP columns
BPDCN712R@meta.data[,which(grepl("MTAP.rearr.\\d", colnames(BPDCN712R@meta.data)))] <- NULL

NormalRelapse <- merge(BM, BPDCN712R)
NormalRelapse$UMAP_1 <- NormalRelapse_old@reductions$umap@cell.embeddings[,1]
NormalRelapse$UMAP_2 <- NormalRelapse_old@reductions$umap@cell.embeddings[,2]
saveRDS(NormalRelapse, file = paste0(dir.ch, "NormalRelapse_Seurat_Genotyped.rds"))


#=============================#
#### Number of transcripts ####
#=============================#

transcripts.tib <- tibble(Mut = unique(FilteredCells_files.tib$Mut), BPDCN628_wt = 0, BPDCN628_mut = 0,
                          BPDCN712_wt = 0, BPDCN712_mut = 0, BPDCN712R_wt = 0, BPDCN712R_mut = 0)

for (n in 1:nrow(FilteredCells_files.tib)) {
    #n <- 1
    UMIs.df <- read.table(FilteredCells_files.tib[n,]$FilteredCells, header = T)
    
    # Filter for unique high-quality cell barcodes
    seu_bcs <- cutf(colnames(get(FilteredCells_files.tib$Sample[n])), d = "-")
    seu_bcs <- setdiff(seu_bcs, seu_bcs[duplicated(seu_bcs)])
    UMIs.df <- UMIs.df[UMIs.df$BC %in% seu_bcs,]
    
    mut <- FilteredCells_files.tib[n,]$Mut
    wt_col <- str_c(FilteredCells_files.tib[n,]$Sample, "_wt")
    mut_col <- str_c(FilteredCells_files.tib[n,]$Sample, "_mut")
    
    transcripts.tib[match(mut, transcripts.tib$Mut), match(wt_col, colnames(transcripts.tib))] <- sum(UMIs.df$wtUMIs)
    transcripts.tib[match(mut, transcripts.tib$Mut), match(mut_col, colnames(transcripts.tib))] <- sum(UMIs.df$mutUMIs)
}

# Save
write_tsv(transcripts.tib, file = str_c(dir.ch, "TranscriptNumbers.txt"))

