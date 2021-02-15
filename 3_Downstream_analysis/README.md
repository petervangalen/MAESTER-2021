The scripts in this folder were used for downstream analysis of the combined scRNA-seq, MAESTER and GoT data.

#### 1. Cell line mixing clustering and cleanup
The same cells, a mixture of K562 and BT142, were analyzed using Seq-Well S^3 and 10X 3' v3 scRNA-seq protocols. These scripts were used to cluster and remove contaminated cells with decontX.\
200915_TenX_CellLineMix_decontX.R\
201101_SW_CellLineMix_decontX.R

#### 2. Cell line mixing variants
These scripts were used to classify cells as either K562 or BT142 based on mitochondrial variants, and to compare this to classification based on RNA-seq.\
201101_SW_CellLineMix_variants.R\
201129_TenX_CellLineMix_variants.R

These scripts were used to identify clonal structure in K562 cells of both cell line mixing experiments.\
201119_SW_K562_clones.R\
201203_TenX_K562_clones.R

#### 3. Patient 10 Diagnosis
These scripts were used to visualize UMAPs, identify informative variants, generate a heatmap of variant allele frequencies, and assess cell type proportions in cells from the uninvolved bone marrow from Patient 10 when he was diagnosed with BPDCN.\
210123_BPDCN712_UMAP.R\
210124_Variants_Of_Interest.R\
210204_LineageBias.R

#### 4. TET2 mutations
This script was used to quantify transcripts with TET2 mutations, that were detected by Genotyping-of-Transcriptomes, in clones that were detected by MAESTER.\
210201_Heatmap.R

#### 5. Cytotoxic T Lymphocyte correlation
This script was used to correlate the transcriptomes of cytotoxic T lymphocytes (CTLs) within (intraclonal) and between (interclonal) clones that were detected by MAESTER.\
210201_CTL_correlation.R

#### 6. Trajectory analysis
This script was used to analyze the myeloid differentiation trajectory including the assignment of pseudotime values and assessment of the density of cells along the trajectory.\
210209_Slingshot.R






