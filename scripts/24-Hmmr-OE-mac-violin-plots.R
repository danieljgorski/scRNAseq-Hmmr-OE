# Violin plots of differentially expressed genes in macrophages

# Load libraries
library(Seurat) #>=4.0.1
library(ggplot2)
library(dplyr)

# Load common downregulated dge data
dge <- read.csv(file = "results/mac-differential-gene-expression/mac_dge_no_threshold.csv")
df <- read.csv(file = "results/mac-common-dge/mac-common-down-dge.csv")

# Load mac subset object for background genes
load("results/objects/mac_annotated.Rdata")

for (i in df$gene) {
  sig <- dge %>%
    filter(gene == i) %>%
    filter(abs(avg_log2FC) > 0.25) %>%
    filter(p_val_adj < 0.01)
  v <- VlnPlot(mac,
               i,
               split.by = "genotype",
               pt.size = 0,
               idents = sig$cluster,
               cols = c("#D4D4D4", "#008000"))
  pdf(file = paste0("results/mac-violin-plots/VlnPlot_", i, ".pdf"),
      useDingbats = F)
  print(v)
  dev.off()
}
