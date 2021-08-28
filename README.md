# MAESTER-2021

## Outline
This repository contains a collection of scripts for the analysis of the MAESTER paper. The following figure (also Supplemental Figure 5 in the paper) shows the general workflow.

![outline](Auxiliary_files/Figure_S5_pipelines.png)



## Pre-processing of sequencing reads and quality controls

Filter for cell barcodes (CBs) and generate fastq files with CB and unique molecular identifiers (UMIs) from the Read 1 fastq in the read ID of the Read 2 fastq:
[Assemble_fastq.R](Pre-processing/Assemble_fastq.R)

After alignment, take the bam file and add the CB and UMI from the read ID as bam tags:
[Tag_CB_UMI.sh](Pre-processing/Tag_CB_UMI.sh)

Process IronThrone-GoT summary tables by doing additional QC and generating tables of wild-type/mutant cells:
[GoT_QC.R](Pre-processing/GoT_QC.R)



## Downstream analyses of combined scRNA-seq, MAESTER and GoT data.

#### 1. Visualize mitochondrial genome coverage
Take maegatk output and plot coverage along the mitochondrial genome

#### 2. Cell line mixing: clustering and cleanup
The same cells, a mixture of K562 and BT142, were analyzed using Seq-Well S^3 and 10X 3' v3 scRNA-seq protocols. These scripts were used to cluster and remove contaminated cells with decontX:

#### 3. Cell line mixing: variants
These scripts in this folder were used to classify cells as either K562 or BT142 based on mitochondrial variants, to compare this to classification based on RNA-seq, and to identify clonal structure in K562 cells of both cell line mixing experiments.

#### 4. Clonal hematopoiesis sample
These scripts were used to analyze cells from the clonal hematopoiesis bone marrow aspirate. This entails UMAP visualization, identification of informative variants, generating a heatmap of variant allele frequencies, and assessing cell type proportions (including the radar plot).

#### 5. TCR-Seq
This folder includes the analyses of the T-cell Receptor Enrichment to linK clonotypes by sequencing protocol. Duncan performed the initial processing using Jupyter, followed by downstream analysis to link mtDNA variants to TCR clonotypes.

#### 6. TET2 mutations
This folder includes analyses to quantify transcripts with TET2 mutations, that were detected by Genotyping of Transcriptomes, in clones that were detected by MAESTER.

#### 7. Trajectory analysis
Slingshot was used to analyze the myeloid differentiation trajectory including the assignment of pseudotime values and assessment of the density of cells along the trajectory.

#### 8. Glioblastoma
We performed Seq-Well scRNA-seq and MAESTER on two regions of a glioblastoma tumor and PBMCs from the same patient. The scripts ultimately leading to the identification of a variant that informs the origin of tumor-infiltrating myeloid cells are here.

#### 9. K562 bulk ATAC-seq
We performed bulk ATAC-sequencing to validate variants that were detected by MAESTER. Here we correlate the VAFs between both methods.




