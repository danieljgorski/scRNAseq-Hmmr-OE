# Volcano plots of differential gene expression analysis of macrophages

# Load libraries
library(Seurat) #>=4.0.1
library(ggplot2)
library(dplyr)
library(ggrepel)
source("scripts/etc/VolcanoPlot.R")

# Load dge data
mac_dge_no_threshold <-
  read.csv(file = "results/mac-differential-gene-expression/mac_dge_no_threshold.csv")

# Loop through clusters and output VolcanoPlots
for (i in unique(mac_dge_no_threshold$cluster)) {
  pdf(
    file = paste0(
      "results/mac-volcano-plots/VolcanoPlot_",
      i,
      ".pdf"
    ),
    height = 6,
    width = 8,
    useDingbats = F
  )
  VolcanoPlot(
    df = mac_dge_no_threshold,
    identity = i,
    title = paste0("Differential gene expression: ", i, " - OE/WT")
  )
  dev.off()
}
