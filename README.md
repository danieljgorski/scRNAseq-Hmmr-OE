# scRNAseq-Hmmr-OE
This repository contains single-cell RNA sequencing analysis of cardiac interstitial cells from WT and Hyaluronan Mediated Motility Receptor(*Hmmr*) over-expressing (OE) mice, 3 days (n = 3, 4) and 7 days (n = 3, 4) after acute myocardial infarction. Data were processed with the Seurat toolkit, using SCTransform normalization and reference-based integration with reciprocal PCA. WT samples were used as reference, OE samples as query.

## Sequencing data
Sequencing data, including fastq files and count matrices will be available upon publication or request.

## Analysis
To recreate the full analysis you can follow the steps below. If you would like to process the data with your own custom workflow, a final list of cells (barcodes + metadata) after preprocessing, doublet removal and low-quality cluster removal can be found in:

* `data/basic_annotation.csv` (all interstitial cells)
* `data/mac_annotation.csv` (re-clustered macrophage subset)

### Libraries
The following scRNA-seq-specific libraries were used:

* [Seurat](https://satijalab.org/seurat/index.html) v4.0.1
* [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) v2.0.3
* [ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap) v2.13.1
* [miloR](https://marionilab.github.io/miloR/) v.1.2.0
* [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) v1.16.0
* [scater](https://bioconductor.org/packages/release/bioc/html/scater.html) v1.22.0
* [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) v3.36.0
* [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) v4.2.1
* [enrichplot](https://bioconductor.org/packages/release/bioc/html/enrichplot.html) v1.14.1
* [org.Mm.eg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html) v3.14.0
* [yulab.utils](https://cran.r-project.org/package=yulab.utils) v0.0.4
* [condiments](https://www.bioconductor.org/packages/release/bioc/html/condiments.html) v1.2.0
* [tradeSeq](https://bioconductor.org/packages/release/bioc/html/tradeSeq.html) v1.8.0
* [slingshot](https://bioconductor.org/packages/release/bioc/html/slingshot.html) v2.2.1

### Instructions
To reproduce the analysis, clone this repository and place the filtered feature count matrices inside the `data` folder. It should then contain the following:
```
scRNAseq-Hmmr-OE/data/
    day_3_filtered_feature_bc_matrix
    day_7_filtered_feature_bc_matrix
    basic_annotation.csv
    canonical_markers.csv
    cell_migration_GO_0016477.csv
    leukocyte_chemotaxis_GO_0030595.csv
    mac_annotation.csv
    seurat_cell_cycle.csv  
```

`0-Hmmr-OE-complete-analysis.R` will create all necessary directories and run the full analysis in the appropriate order. Each analysis step can also be run individually for better interactivity, starting from `1-Hmmr-OE-preprocessing.R`. The estimated total memory needed to hold and process the resulting Seurat objects is ~120GB.

## Examples
<p align="center">
  <img src="/examples/DimPlot_basic_annotation.png" width="1000">
</p>

<p align="center">
  <img src="/examples/Heatmap.png" width="1000">
</p>


## To-do
* Add abstract at submission
* Add full author list at submission
* Add differential expression + mechanism schematic after publication
