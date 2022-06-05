## scRNAseq-Hmmr-OE

### Overview
Single-cell RNA sequencing of 137,166 cardiac interstitial cells from WT and Hyaluronan Mediated Motility Receptor(*Hmmr*) over-expressing (OE) mice, 3 (n = 3, 4) and 7 days (n = 3, 4) after acute myocardial infarction. Data were processed with the Seurat toolkit, using SCTransform normalization and reference-based integration with reciprocal PCA. WT samples were used as reference, OE samples as query.

### Abstract

### Sequencing data
Sequencing data, including fastq files and count matrices will be available upon publication.

### Libraries
Seurat v4.0.1
DoubletFinder v2.0.3
dplyr
ggplot2
library(patchwork)
library(readr)
library(ggrepel)
library(ComplexHeatmap) v2.13.1
library(stringr)
library(scProportionTest)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(edgeR)
library(yulab.utils)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(knitr)
library(kableExtra)
library(forcats)

### Analysis
`0-Hmmr-OE-complete-analysis.R` will create all necessary directories and run the full analysis in the appropriate order. Each analysis step can also be run individually for better interactivity. The estimated total memory needed to hold and process the resulting Seurat objects is ~120GB. Scripts 1-6 perform preprocessing, doublet detection, low quality cluster removal and reference-based integration. If you would like to skip these and cluster the data on your own, using only the final cells, their barcodes and necessary metadata are stored in `data/basic_annotation.csv` (all interstitial cells) and `data/mac_annotation.csv`(macrophage subset re-clustering).

### Examples
<p align="center">
  <img src="/eg/DimPlot_basic_annotation.png" width="1000">
</p>

<p align="center">
  <img src="/eg/Heatmap.png" width="1000">
</p>
