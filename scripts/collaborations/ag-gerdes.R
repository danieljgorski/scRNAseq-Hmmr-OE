# Collboration with AG Gerdes - Alexander Lang - CD40, CD40 ligand

# Load libraries
library(Seurat) #>=4.0.1
library(ggplot2)
library(ggrepel)
source("scripts/etc/colors.R")

# Load object, subset for WT cells
load("results/objects/obj_annotated.Rdata")
obj <- subset(obj, subset = genotype == "WT")

# GOIs
GOI <- c("Cd40", "Cd40lg", "Ccl5", "Ccr7")

# Annotated Dimplot of WT cells
p <- DimPlot(obj,
             reduction = "umap",
             pt.size = .5,
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
jpeg(file = "results/collaborations/ag-gerdes/DimPlot.jpg",
     height = 2200,
     width = 2500,
     res = 300,
     units = "px")
print(q)
dev.off()

# FeaturePlots
for (i in GOI) {
  p <- FeaturePlot(obj, 
                   features = i, 
                   pt.size = .5,
                   cols = c("lightgrey", "red"),
                   reduction = "umap",
                   raster = F) +
    ggtitle(i) +
    theme(legend.position = c(0.9, 0.15),
          plot.title = element_text(face = "italic", size = 35),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  jpeg(file = paste0("results/collaborations/ag-gerdes/FeaturePlot_", i,".jpg"),
       height = 2200,
       width = 2500,
       res = 300,
       units = "px")
  print(p)
  dev.off()
}

# FeaturePlots by day
for (i in GOI) {
  p <- FeaturePlot(obj, 
                   features = i, 
                   pt.size = .5,
                   cols = c("lightgrey", "red"),
                   reduction = "umap",
                   raster = F,
                   split_)
  jpeg(file = paste0("results/collaborations/ag-gerdes/FeaturePlot_", i,".jpg"),
       height = 2200,
       width = 2500,
       res = 300,
       units = "px")
  print(p)
  dev.off()
}

# VlnPlot
for (i in GOI) {
  p <- VlnPlot(obj,
               features = i,
               pt.size = 0.65) + 
    NoLegend() + 
    ylab(paste0(i, " expression")) +
    theme(plot.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14)) +
    RotatedAxis()
  jpeg(file = paste0("results/collaborations/ag-gerdes/VlnPlot_", i, ".jpg"),
       height = 1100,
       width = 4000,
       res = 300,
       units = "px")
  print(p)
  dev.off()
}

# VlnPlot by day
for (i in GOI) {
  p <- VlnPlot(obj,
               features = i,
               pt.size = 0.65,
               split.by = "timepoint") + 
    ylab(paste0(i, " expression")) +
    theme(plot.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14)) +
    RotatedAxis()
  jpeg(file = paste0("results/collaborations/ag-gerdes/VlnPlot_", i, "_by_day.jpg"),
       height = 1100,
       width = 8000,
       res = 300,
       units = "px")
  print(p)
  dev.off()
}

# VlnPlot by day
for (i in GOI) {
  p <- VlnPlot(obj,
               features = i,
               pt.size = 0.65,
               split.by = "timepoint",
               cols = c("#ED553B", "#20639B")) + 
    ylab(paste0(i, " expression")) +
    theme(plot.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14)) +
    RotatedAxis()
  jpeg(file = paste0("results/collaborations/ag-gerdes/VlnPlot_", i, "_by_day.jpg"),
       height = 1100,
       width = 7000,
       res = 300,
       units = "px")
  print(p)
  dev.off()
}

# DotPlot
jpeg(file = paste0("results/collaborations/ag-gerdes/DotPlot.jpg"),
     height = 800,
     width = 6000,
     res = 300,
     units = "px")
p <- DotPlot(obj,
        features = GOI,
        split.by = "timepoint",
        group.by = rev("basic_annotation"),
        cols = c("#20639B", "#ED553B")) +
  coord_flip() +
  RotatedAxis() +
  scale_y_discrete(limits = rev) +
  scale_x_discrete(limits = rev)
print(p)
dev.off()

