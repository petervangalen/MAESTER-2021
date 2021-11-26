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
Take maegatk output and plot coverage along the mitochondrial genome.

#### 2. Cell line mixing: clustering and cleanup
The same cells, a mixture of K562 and BT142, were analyzed using Seq-Well S^3 and 10x 3' v3 scRNA-seq protocols. These scripts were used to cluster and remove contaminated cells with decontX.

#### 3. Cell line mixing: variants
The scripts in this folder were used to classify cells as either K562 or BT142 based on mitochondrial variants, to compare this to classification based on RNA-seq, and to identify clonal structure in K562 cells of both cell line mixing experiments.

#### 4. Clonal hematopoiesis sample
These scripts were used to analyze cells from the clonal hematopoiesis bone marrow aspirate. This entails UMAP visualization, selection of informative variants, generating heatmaps of variant allele frequencies, and assessing cell type proportions (including the radar plot).

#### 5. TCR-Seq
This folder includes the analyses of the T-cell Receptor Enrichment to linK clonotypes by sequencing (TREK-seq) protocol. Duncan Morgan performed initial processing using Jupyter, followed by downstream analyses to link mtDNA variants to TCR clonotypes.

#### 6. TET2 mutations
This folder includes analyses to intersect TET2 mutations (detected by Genotyping of Transcriptomes) with mtDNA variants (detected by MAESTER).

#### 7. Trajectory analysis
Slingshot was used to analyze HSC differentiation trajectories, including the assignment of pseudotime values and density plots of cells along the trajectories.

#### 8. Glioblastoma
We performed Seq-Well scRNA-seq and MAESTER on peripheral blood mononuclear cells and glioblastoma tumor regions from a patient at the time of surgery. This includes selection of informative variants, the identification of a tumor-specific deletion, and the identification of a variant that suggests a blood origin of some tumor-infiltrating myeloid cells.

#### 9. K562 bulk ATAC-seq
We performed bulk ATAC-sequencing to validate variants that were detected by MAESTER. Here we correlate the VAFs between both methods.




