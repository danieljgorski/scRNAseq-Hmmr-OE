# Cluster markers, heatmap and average expression of annotated data

# Load libraries
library(Seurat) #>=4.0.1
library(ggplot2)
library(dplyr)
library(ggrepel)
source("scripts/colors.R")

# Load object
load("results/objects/obj_annotated.Rdata")

# Find cluster markers
Idents(obj) <- "basic_annotation"
cluster_markers <- FindAllMarkers(obj,
                                  assay = "RNA",
                                  logfc.threshold = 0.5,
                                  only.pos = T,
                                  return.thresh = 0.001,
                                  densify = T,
                                  verbose = T)
write.csv(cluster_markers,
          file = "results/cluster-markers/cluster_markers.csv",
          row.names = F)

# Heatmap
top20 <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
Idents(obj) <- "basic_annotation"
p <- DoHeatmap(subset(obj, downsample = 500), 
               features = top20$gene,
               label = F,
               assay = "RNA",
               group.colors = colors,
               group.by = "basic_annotation") +
  scale_fill_gradient2(low = "#2664ad", mid = "white", high = "#eb2a0e") +
  ylab("Marker genes") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_text(size = 18),
        legend.text= element_text(size= 12),
        legend.title = element_text(size = 18),
        legend.justification = "top") +
  scale_color_manual(values = colors) +
  guides(color = guide_legend(title = "Identity",
                              override.aes = list(shape=16,
                                                  size=4.25,
                                                  alpha=1)))
pdf(file = "results/cluster-markers/Heatmap_top_20.pdf",
    height = 9,
    width = 7,
    useDingbats = F)
print(p)
dev.off()

# Return averaged expression values for each identity class, convert to df
avg_exp <- AverageExpression(obj, assays = "RNA", group.by = "basic_annotation")
avg_exp <- avg_exp$RNA
avg_exp <- as.data.frame(avg_exp)

# Loop through each cluster, find top 20 highest expressed genes, save
avg_exp_list <- list()
for (i in colnames(avg_exp)) {
  top_20 <- avg_exp %>% select(i) %>% top_n(20) 
  top_20$gene <- rownames(top_20)
  top_20$cluster <- i
  colnames(top_20) <- c("expression", "gene", "cluster")
  top_20 <- top_20[, c("cluster", "gene", "expression")]
  avg_exp_list[[i]] <- top_20
}
avg_exp_top_20 <- do.call(rbind, avg_exp_list)
write.csv(avg_exp_top_20,
          file = "results/cluster-markers/avg_exp_top_20.csv",
          row.names = F)
