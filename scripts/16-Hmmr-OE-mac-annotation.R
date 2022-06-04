# Macrophage annotation

# Load libraries
library(Seurat) #v4.0.1
library(dplyr)
library(ComplexHeatmap)
library(ape)
library(ggplot2)
library(patchwork)
library(yulab.utils)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
source("scripts/etc/ORA.R")

# Load full object
load("results/objects/mac_integrated_clean.Rdata")

# Evaluation of un-annotated seurat clusters

# Dimplots broken down by time point----
p1 <- DimPlot(mac, label = T) + NoLegend()
p2 <- DimPlot(mac, label = T, split.by = "timepoint")
pdf(file = "results/mac-annotation/seurat_cluster_dimplots.pdf",
    height = 6,
    width = 12,
    useDingbats = F)
p1 + p2
dev.off()

# Calculate cluster markers----
seurat_cluster_markers_mac <- FindAllMarkers(mac,
                                         assay = "RNA",
                                         only.pos = T,
                                         densify = T,
                                         verbose = T)
write.csv(seurat_cluster_markers_mac,
          file = "results/mac-annotation/seurat_cluster_markers_mac.csv",
          row.names = F)
top5 <- seurat_cluster_markers_mac %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# Heatmap----

# Filter top 20 marker genes, and top 1 for labeling
top20 <- seurat_cluster_markers_mac %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
top20 <- top20[!duplicated(top20$gene), ]
top2 <- top20 %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)
top2 <- top2[!duplicated(top2$gene), ]

# Filter out duplicates
top20 <- top20[!duplicated(top20$gene), ]
top2 <- top2[!duplicated(top2$gene), ]

# Down sample extract markers from scaled data slot
mac_down <- subset(mac, downsample = 100)
mat <- mac_down[["RNA"]]@scale.data[top20$gene, ] %>% as.matrix()

# Find scaled data range and create custom colors
quantile(mat, c(0.05, 0.95))
cols <- circlize::colorRamp2(c(-1.5, 0, 2), c("#2664ad", "white", "#eb2a0e"))

# Build annotations for column and row splitting
cluster_anno <- mac_down@meta.data$seurat_clusters
right <- rowAnnotation(foo = anno_mark(which(rownames(mat) %in%
                                               top2$gene),
                                       labels = top2$gene,
                                       labels_gp = gpar(fontsize = 10,
                                                        fontface = "italic")))

# Plot
p <- Heatmap(mat,
             name = "Scaled expression",
             col = cols,
             column_split = cluster_anno,
             column_gap = unit(0.6, "mm"),
             cluster_columns = F,
             show_column_names = F,
             column_title_gp = gpar(fontsize = 14),
             column_title_rot = 90,
             right_annotation = right,
             row_split = factor(top20$cluster, levels = levels(cluster_anno)),
             row_gap = unit(0.6, "mm"),
             cluster_rows = F,
             show_row_names = F,
             row_title = NULL,
             use_raster = T,
             raster_quality = 4,
             heatmap_width = unit(24, "cm"),
             heatmap_height = unit(27, "cm"),
             heatmap_legend_param = list(direction = "horizontal",
                                         title_gp = gpar(fontize = 2),
                                         title_position = "lefttop"))
pdf(file = "results/mac-annotation/seurat_cluster_heatmap.pdf",
    width = 10,
    height = 12,
    useDingbats = F)
draw(p, heatmap_legend_side = "bottom")
dev.off()

# Canonical marker gene expression----

# Macrophage genes
mac_genes <- c("Ly6c2",
               "Ccr2",
               "Cx3cr1",
               "Il1b",
               "Il6",
               "Tnf",
               "Il10",
               "Tgfb1",
               "Adgre1",
               "H2-Aa",
               "Itgam",
               "Lgals3",
               "Fcgr1",
               "Mrc1",
               "Maf",
               "Trem2",
               "Mertk",
               "C1qa",
               "C1qb",
               "C1qc",
               "Igf1",
               "Pdgfb",
               "Pdgfc",
               "Timd4",
               "Lyve1",
               "S100a9",
               "Csf3r",
               "Ly6g",
               "Siglecf",
               "Retnlg",
               "G0s2",
               "Lcn2",
               "Wfdc21",
               "Ifit3",
               "Ifit1",
               "Cxcl10",
               "Rsad2",
               "Ear2",
               "Pglyrp1",
               "Eno3",
               "Itgal",
               "Ace",
               "Saa3",
               "Ltc4s",
               "Ccl24",
               "Ednrb",
               "Mki67",
               "Ccnb2")

# Loop through markers and generate Feature and VlnPlots of expression
for (i in mac_genes) {
  p1 <- FeaturePlot(mac, features = i, label = T, raster = F)
  p2 <- VlnPlot(mac, features = i, pt.size = 0.01, sort = T) + NoLegend() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  pdf(file = paste0("results/mac-annotation/Feat_VlnPlot_", i, ".pdf"),
      width = 16,
      height = 5,
      useDingbats = F)
  print(p1 + p2 + plot_layout(ncol = 2, widths = c(1, 2)))
  dev.off()
}


# Gene ontology over-representation analysis of marker genes----

# Loop through clusters, extract markers and input into GOBP ora
for (i in unique(seurat_cluster_markers_mac$cluster)) {
  df <- seurat_cluster_markers_mac %>% dplyr::filter(cluster == i)
  pdf(file = paste0("results/mac-annotation/Cluster_", i, "_ora.pdf"),
                    useDingbats = F)
  ORA(df$gene, rownames(mac), paste0("Cluster-", i, ": cluster markers"))
  dev.off()
}

###############################################################################
# Annotation notes
###############################################################################
# 8-(Mac-TR, tissue resident macrophages, d3 + d7, Ly6c2-, Ccr2-, Adgre1+
#   Cx3cr1+ Timd4+, Lyve1+, [Sub-population is Timd4- Lyve1- H2-Aa-high,
#   similar to Farbehi et al. 2019])
#
# 6-(Mono, bone marrow derived monocytes, d3, Ly6c2+, Ccr2-high, Adgre1-low,
#   H2-Aa-low, Fcgr1+, Il1b-high, [Low expression of typical macrophage
#   functional markers, phagocytosis and complement, i.e. Mrc1-low, Trem2-low,
#   Mertk-low, C1q(a-c)-low])
#
# 2-(Mac-Ifn, interferon stimulated macrophages, Ly6c2+, Ccr2-high, Adgre1+,
#    Ifit1+, Ifit3+, Cxcl10+, Rsad2+)
#
# 0-(Mac-M1-1, monocyte derived pro-inflammatory classical M1 macrophages,
#   largely d3, some d7, Ly6c2-, Ccr2+, Adgre1+, H2-Aa-low, Cx3cr1-low, Mrc1+,
#   Trem2+, Mertk+, C1q(a-c)-high)
#
# 3-(Mac-M1-2, monocyte derived pro-inflammatory classical M1 macrophages,
#   largely d3, Ly6c2-, Ccr2+, Adgre1+, H2-Aa-low, Cx3cr1-, Mrc1+, Trem2+,
#   Mertk+, C1q(a-c)-high)
#
# 5-(Mac-M1-3, monocyte derived pro-inflammatory classical M1 macrophages,
#   largely d3, Ly6c2-, Ccr2+, Adgre1+, H2-Aa-low, Cx3cr1-, Mrc1+, Trem2+,
#   Mertk+, C1q(a-c)-high)
#
# 1-(Mac-M2-1, monocyte derived reparative non-classical M2 macrophages,
#   d3 and d7, Ly6c2-, Ccr2-high, Adgre1+, H2-Aa-high, Cx3cr1+)
#
# 11-(Mac-M2-2, monocyte derived reparative non-classical M2 macrophages,
#   largely d3, Ly6c2-low, Ccr2-high, Adgre1+, H2-Aa-high, Cx3cr1+, [May be
#   directly differentiated from monocyte pool? Ly6c2-low (not negative),
#   Il1b-high, low expression of typical macrophage functional markers,
#   phagocytosis and complement, i.e. Mrc1-low, Trem2-low, Mertk-low,
#   C1q(a-c)-low])
#
# 13-(Mac-M2-3, monocyte derived reparative non-classical M2 macrophages,
#   d3 and d7, Ly6c2- Ccr2+, Adgre1-low, H2-Aa-high)
#
# 7-(Mac-Cyc-1, cycling macrophages, Ly6c2- Ccr2+ Adgre1+ Cx3cr1+ Mki67+ Ccnb2+)
#
# 10-(Mac-Cyc-2, cycling macrophages, Ly6c2- Ccr2+ Adgre1+ Cx3cr1+ Mki67+)
#
# 4-(Mac-12, Lgals3+, Fcgr1+, Trem2+, Ly6c2-, Ccr2-, Adgre1-, Mrc1-low,
#    Mertk-low)
#
# 9-(Mac-13, similar to MAC8 from Farbehi et al. 2019, Ly6c2-, Ccr2+, Adgre1+,
#   Saa3+, Ltc4s+, Ccl24+, Ednrb+)
# 12-(Mac-14, Ly6c2-, Ccr2+, Adgre1+, Mrc1-low, Mertk-low)
#
# 14-(Mac-15, similar to MAC7 from Farbehi et al. 2019, Ly6c2-, Ccr2+, Adgre1+,
#   Cx3cr1+, C1qa-c-Low, Mrc1-low, Trem2-low, Mertk-low, Il1b-high, Ear2+,
#   Pglyrp1+, Eno3+, Itgal+, Ace+)
#
# 15-(Mac-16, Ly6c2-, Ccr2+, Adgre1+, H2-Aa+)
###############################################################################

# Rename Idents to annotations
mac <- RenameIdents(mac,
                    "0" = "Mac-M1-1",
                    "1" = "Mac-M2-1",
                    "2" = "Mac-IFN",
                    "3" = "Mac-M1-2",
                    "4" = "Mac-12",
                    "5" = "Mac-M1-3",
                    "6" = "Mono",
                    "7" = "Mac-Cyc-1",
                    "8" = "Mac-TR",
                    "9" = "Mac-13",
                    "10" = "Mac-Cyc-2",
                    "11" = "Mac-M2-2",
                    "12" = "Mac-14",
                    "13" = "Mac-M2-3",
                    "14" = "Mac-15",
                    "15" = "Mac-16")

# Store renamed idents as a new meta data column
mac@meta.data$mac_annotation <- Idents(mac)

# Refactor annotation levels
source("scripts/etc/maclevels.R")
mac@meta.data$mac_annotation <- factor(mac@meta.data$mac_annotation,
                                         levels = maclevels)
DimPlot(mac,
        group.by = "mac_annotation",
        label = T,
        repel = T)

# Set Idents as re-factored mac_annotation identities
Idents(mac) <- "mac_annotation"

# Save object with basic annotations
save(mac, file = "results/objects/mac_annotated.Rdata")

# Save annotation, barcodes and UMAP embeddings etc. for consistency in
# external usage, de-comment to overwrite
# barcodes <- rownames(mac@meta.data)
# annotation <- mac@meta.data$mac_annotation
# genotype <- mac@meta.data$genotype
# timepoint <- mac@meta.data$timepoint
# sample <- mac@meta.data$sample
# UMAP_1 <- Embeddings(mac[["umap"]])[,1]
# UMAP_2 <- Embeddings(mac[["umap"]])[,2]
# mac_annotation <- data.frame(barcodes,
#                                annotation,
#                                genotype,
#                                timepoint,
#                                sample,
#                                UMAP_1,
#                                UMAP_2)
# write.csv(mac_annotation, file = "data/mac_annotation.csv", row.names = F)
