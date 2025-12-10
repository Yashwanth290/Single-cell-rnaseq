Single-Cell RNA-seq Analysis of SARS-CoV-2 Infection

it contains code and workflows for single-cell RNA-seq (scRNA-seq) analysis aimed at characterizing cellular and transcriptional responses to SARS-CoV-2 infection.
The analysis includes quality control, normalization, clustering, differential gene expression (DGE), and pathway enrichment across infected and uninfected cells.


Analysis Overview
1. Data Source

Single-cell transcriptomes from SARS-CoV-2 infected and control samples

2. Preprocessing & Quality Control

Includes:

Import of Cell Ranger matrices or FASTQ processing

Filtering using thresholds such as:

min.cells

min.features

percent.mt (mitochondrial content)

Doublet detection (optional: DoubletFinder or Scrublet)

QC metrics saved to results/figures/qc/

3. Normalization

Supports either:

SCTransform (default for Seurat-based workflow)

LogNormalize

Highly variable gene (HVG) selection

4. Dimensionality Reduction & Clustering

PCA

Nearest neighbor graph

UMAP / t-SNE

Clustering using the Leiden/Louvain algorithm

5. Cell Type Annotation

Manual curation from markers

Automated annotation tools (optional):

SingleR, scCATCH, CellTypist, Azimuth

6. Differential Expression Analysis

DGE performed:

Between infected vs non-infected cells

Between cell types

For cluster markers

Methods include:

Seurat::FindMarkers() (Wilcoxon, MAST, DESeq2, logistic regression)

Exported to results/dge_tables/

7. Pathway Analysis

Gene set enrichment using:

GO: Biological Process

KEGG pathways

Reactome signaling

MSigDB Hallmark pathways

Requirements
R Version

≥ 4.2 recommended

Essential Packages

Seurat

dplyr

patchwork

harmony 

DoubletFinder 

clusterProfiler

msigdbr

ggplot2
