# VolcanoPlots of differential gene expression analysis

# Load libraries
library(Seurat) #>=4.0.1
library(ggplot2)
library(dplyr)
library(ggrepel) 
source("scripts/VolcanoPlot.R")

# Load dge data
dge_no_threshold <- 
  read.csv(file = "results/differential-gene-expression/dge_no_threshold.csv")

# Loop through clusters and output VolcanoPlots
for (i in unique(dge_no_threshold$cluster)) {
  pdf(
    file = paste0(
      "results/volcano-plots/VolcanoPlot_",
      i,
      ".pdf"
    ),
    height = 6,
    width = 8,
    useDingbats = F
  )
  VolcanoPlot(
    df = dge_no_threshold,
    identity = i,
    title = paste0("Differential gene expression: ", i, " - OE/WT")
  )
  dev.off()
}