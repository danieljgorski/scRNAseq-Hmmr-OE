# Over-representation analysis of differentially expressed genes in macrophages

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
source("scripts/ORA_up.R")
source("scripts/ORA_down.R")

# Load dge data
mac_dge_no_threshold <-
  read.csv(file = "results/mac-differential-gene-expression/mac_dge_no_threshold.csv")

# Load object for background genes
load("results/objects/mac_annotated.Rdata")
detected <- rownames(mac)
remove(mac)

# Loop through clusters and export upregulated ORA analyses
for (i in unique(mac_dge_no_threshold$cluster)) {
  pdf(file = paste0("results/mac-ora/ORA_up_", i, ".pdf"),
      height = 4,
      width = 7,
      useDingbats = F)
  ORA_up(mac_dge_no_threshold, i, detected, "Upregulated in OE")
  dev.off()
}

# Loop through clusters and export downregulated ORA analyses
for (i in unique(mac_dge_no_threshold$cluster)) {
  pdf(file = paste0("results/mac-ora/ORA_down_", i, ".pdf"),
      height = 4,
      width = 7,
      useDingbats = F)
  ORA_down(mac_dge_no_threshold, i, detected, "Downregulated in OE")
  dev.off()
}
