# Cluster markers, heatmap and average expression of annotated macrophages

# Load libraries
library(Seurat) #>=4.0.1
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ComplexHeatmap)
source("scripts/etc/colors.R")
source("scripts/etc/maclevels.R")

# Load object
load("results/objects/mac_annotated.Rdata")

# Find cluster markers
mac_cluster_markers <- FindAllMarkers(mac,
                                  assay = "RNA",
                                  logfc.threshold = 0.5,
                                  only.pos = T,
                                  return.thresh = 0.001,
                                  densify = T,
                                  verbose = T)
write.csv(mac_cluster_markers,
          file = "results/mac-cluster-markers/mac_cluster_markers.csv",
          row.names = F)

# Heatmap

# Filter top 20 marker genes, and top 1 for labeling
top20 <- mac_cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
top20 <- top20[!duplicated(top20$gene),]
top3 <- top20 %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC)
top3 <- top3[!duplicated(top3$gene),]

# Filter out duplicates
top20 <- top20[!duplicated(top20$gene),]
top3 <- top3[!duplicated(top3$gene),]

# Extract markers from scaled data slot
mat <- mac[["RNA"]]@scale.data[top20$gene, ] %>% as.matrix()

# Find scaled data range and create custom colors
quantile(mat, c(0.05, 0.95))
cols = circlize::colorRamp2(c(-2, 0, 2), c("#2664ad", "white", "#eb2a0e"))

# Build annotations for column and row splitting
cluster_anno <- mac@meta.data$mac_annotation
cluster_anno <- factor(cluster_anno, levels = maclevels)
top <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = colors[1:16]),
                                          height = unit(4, "mm")))
right <- rowAnnotation(foo = anno_block(gp = gpar(fill = colors[1:16]),
                                        width = unit(2, "mm")),
                       bar = anno_mark(which(rownames(mat) %in%
                                                   top3$gene),
                                           labels = top3$gene,
                                           labels_gp = gpar(fontsize = 10,
                                                            fontface = "italic")))

# Plot heatmap
p <- Heatmap(mat,
        name = "Scaled expression",
        col = cols,
        top_annotation = top,
        column_split = cluster_anno,
        column_gap = unit(0.6, "mm"),
        cluster_columns = F,
        show_column_names = F,
        column_title_gp = gpar(fontsize = 14),
        column_title_rot = 90,
        right_annotation = right,
        row_split = factor(top20$cluster, levels = maclevels),
        row_gap = unit(0.6, "mm"),
        cluster_rows = F,
        show_row_names = F,
        row_title = NULL,
        use_raster = T,
        raster_quality = 4,
        heatmap_width = unit(24, "cm"),
        heatmap_height = unit(24, "cm"),
        heatmap_legend_param = list(direction = "horizontal",
                                    title_gp = gpar(fontize = 2),
                                    title_position = "lefttop"))
pdf(file = "results/mac-cluster-markers/mac_heatmap.pdf",
    width = 10,
    height = 12,
    useDingbats = F)
draw(p, heatmap_legend_side = "bottom")
dev.off()

# Return averaged expression values for each identity class, convert to a df
avg_exp <- AverageExpression(mac, assays = "RNA", group.by = "mac_annotation")
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
mac_avg_exp_top_20 <- do.call(rbind, avg_exp_list)
write.csv(mac_avg_exp_top_20,
          file = "results/mac-cluster-markers/mac_avg_exp_top_20.csv",
          row.names = F)
