# DimPlots of annotated macrophage subset
# Load libraries
library(Seurat) #>=4.0.1
library(ggplot2)
library(ggrepel)
source("scripts/etc/colors.R")
source("scripts/etc/HighlightedDimPlot.R")

# Load object
load("results/objects/mac_annotated.Rdata")

# DimPlot of mac annotation
p <- DimPlot(mac,
             reduction = "umap",
             pt.size = .3,
             raster = F,
             label = F,
             cols = colors,
             group.by = "mac_annotation") +
  xlab("UMAP-1") +
  ylab("UMAP-2") +
  labs(color = "Identity") +
  theme(legend.text = element_text(size = 13),
        legend.title = element_text(size = 16),
        legend.justification = "top",
        legend.key.size = unit(3, "point"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 16),
        plot.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4.25),
                              nrow = 37))
q <- LabelClusters(plot = p,
                   id = "mac_annotation",
                   repel = T,
                   force = 0.5,
                   box = T,
                   fill = alpha("white", 0.45),
                   size = 4,
                   label.r = unit(0.25, "lines"),
                   label.size = NA)
pdf(file = "results/mac-dimplots/DimPlot_mac_annotation.pdf",
    height = 6.5,
    width = 8,
    useDingbats = F)
print(q)
dev.off()

# Highlighted DimPlots of each cluster
for (i in levels(Idents(mac))) {
  pdf(file = paste0("results/mac-dimplots/DimPlot_highlighted_", i, ".pdf"),
      height = 6.5,
      width = 8,
      useDingbats = F)
  HighlightedDimPlot(mac, i)
  dev.off()
}
