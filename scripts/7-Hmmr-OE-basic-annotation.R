# Basic cluster annotation based on canonical cell-type markers

# Load libraries
library(Seurat) #v4.0.1
library(patchwork)
library(ggplot2)

# Load object
load("results/objects/obj.Rdata")

# Read in markers
markers <- read.csv(file = "data/canonical_markers.csv")
markers <- markers$All

# Loop through markers and generate Feature and VlnPlots of expression
for (i in markers) {
  p1 <- FeaturePlot(obj, features = i, label = T, raster = F)
  p2 <- VlnPlot(obj, features = i, pt.size = 0, sort = T) + NoLegend() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  pdf(file = paste0("results/basic-annotation/Feat_VlnPlot_", i, ".pdf"),
      width = 16,
      height = 5,
      useDingbats = F)
  print(p1 + p2 + plot_layout(ncol = 2, widths = c(1,2)))
  dev.off()
}


###############################################################################
# Annotation notes
###############################################################################
# 0 - Cdh5, Pecam1, Kdr
# 1 - Fcgr1, Adgre1, Cd68
# 2 - Fcgr1, Adgre1, Cd68, Cd74, H2-Ab1, H2-Aa
# 3 - Fcgr1, Adgre1, Cd68
# 4 - Cd68 low, S100a8, S100a8
# 5 - Fcgr1, Adgre1, Cd68, Cd74, H2-Ab1, H2-Aa
# 6 - Fcgr1, Adgre1, Cd68, Cd74 low
# 7 - Cdh5, Pecam1, Kdr
# 8 - Cd74, H2-Ab1, H2-Aa, Ms4a1, Cd79a
# 9 - Tcf21, Pdgfra, Cthrc1, Postn
# 10 - Cd68
# 11 - Fcgr1, Adgre1, Cd68, Cd74
# 12 - Cd68 low, S100a8
# 13 - Tcf21, Pdgfra, Postn low
# 14 - Fcgr1, Adgre1, Cd68
# 15 - Fcgr1, Adgre1, Cd68, Mki67, Ccnb2, Cd74 low
# 16 - Fcgr1, Adgre1, Cd68, Cd74 low
# 17 - Cdh5, Pecam1, Kdr
# 18 - Fcgr1, Cd68, Cd74, H2-Ab1, H2-Aa, Cd209a DENDRITIC CELLS
# 19 - Fcgr1, Adgre1, Cd68, Cd74 low
# 20 - Mki67 low, Ccnb2 low, Tcf21, Pdgfra, Cthrc1, Postn, Clu, Wt1
# 21 - Fcgr1, Adgre1, Cd68, Cd74 low
# 22 - 
# 23 - Fcgr1, Adgre1, Cd68, Cd74, H2-Ab1, H2-Aa, Postn low, S100a8
# 24 - Cthrc1 low, Postn low
# 25 - Ccl5
# 26 - Ccl5
# 27 - Cd68, Mki67 low, Cd74, H2-Ab1, H2-Aa
# 28 - Mki67, Ccnb2, Cdh5, Pecam1, Kdr, Wt1
# 29 - Clu
# 30 - , Postn low
# 31 - Cd68 low, Cd74, H2-Ab1, H2-Aa, Ccl5
# 32 - Fcgr1, Adgre1, Cd68, Cd74, H2-Ab1, H2-Aa, Cdh5, Pecam1, Kdr
# 33 - Adgre1 low, Cd68, Mki67 low, Ccnb2 low
# 34 - Cdh5, Pecam1, Kdr, Clu
# 35 - Fcgr1, Adgre1, Cd68, Cd74, H2-Ab1 low, H2-Aa low, S100a8, S100a8
# 36 - Cd74, H2-Ab1, H2-Aa, Cdh5, Pecam1, Kdr, Ms4a1, Cd79a
# 37 - Mki67, Ccnb2, Tcf21, Pdgfra, Postn
# 38 -
