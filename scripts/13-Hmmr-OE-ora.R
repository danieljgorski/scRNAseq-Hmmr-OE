# Over-representation analysis of differentially expressed genes

# Load libraries
library(Seurat) #>=4.0.1
library(ggplot2)
library(dplyr)
library(ggrepel)
library(yulab.utils)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(knitr)
library(kableExtra)
source("scripts/etc/ORA_up.R")
source("scripts/etc/ORA_down.R")

# Load dge data
dge_no_threshold <-
  read.csv(file = "results/differential-gene-expression/dge_no_threshold.csv")

# Load object for background genes
load("results/objects/obj_annotated.Rdata")
detected <- rownames(obj)
remove(obj)

# Loop through clusters and export upregulated ORA analyses
for (i in unique(dge_no_threshold$cluster)) {
  pdf(file = paste0("results/ora/ORA_up_", i, ".pdf"),
      height = 4,
      width = 7,
      useDingbats = F)
  ORA_up(dge_no_threshold, i, detected, "Upregulated in OE")
  dev.off()
}

# Loop through clusters and export downregulated ORA analyses
for (i in unique(dge_no_threshold$cluster)) {
  pdf(file = paste0("results/ora/ORA_down_", i, ".pdf"),
      height = 4,
      width = 7,
      useDingbats = F)
  ORA_down(dge_no_threshold, i, detected, "Downregulated in OE")
  dev.off()
}
