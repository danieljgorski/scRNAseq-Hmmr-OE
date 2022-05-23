# DimPlots of annotated data

# Load libraries
library(Seurat) #>=4.0.1
library(ggplot2)
library(ggrepel)
source("scripts/colors.R")

# Load object
load("results/objects/obj_annotated.Rdata")

# DimPlot of basic annotation
p <- DimPlot(obj, 
             reduction = "umap", 
             pt.size = .3,
             raster = F, 
             label = F,
             cols = colors,
             group.by = "basic_annotation") + 
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
                   id = "basic_annotation", 
                   repel = T,
                   force = 0.5,
                   box = T, 
                   fill = alpha("white", 0.45),
                   size = 4,
                   label.r = unit(0.25, "lines"),
                   label.size = NA)
pdf(file = "results/dimplots/DimPlot_basic_annotation.pdf",
    height = 6.5,
    width = 8,
    useDingbats = F)
print(q)
dev.off()

# Highlighted Dimplots
for (i in levels(Idents(obj))) {
  dd <- DimPlot(obj, 
                reduction = "umap", 
                cells.highlight = WhichCells(obj, idents = i),
                cols.highlight = "#5179A7",
                raster = F,
                pt.size = .3) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank()) +
    NoLegend()
  dd <- LabelClusters(plot = dd, 
                      id = "ident",
                      clusters = i,
                      repel = F, 
                      box = T, 
                      fill = alpha("white", 0.45),
                      size = 4,
                      label.r = unit(0.25, "lines"),
                      label.size = NA)
  pdf(file = paste0("results/dimplots/DimPlot_highlighted_", i, ".pdf"),
      height = 6.5,
      width = 8,
      useDingbats = F)
  print(dd)
  dev.off()
}
