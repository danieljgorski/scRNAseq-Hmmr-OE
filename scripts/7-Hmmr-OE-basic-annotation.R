# Basic cluster annotation based on canonical cell-type markers

# Load libraries
library(Seurat) #v4.0.1
library(patchwork)
library(ggplot2)
library(dplyr)
library(ggrepel)

# Load object
load("results/objects/obj_integrated_clean.Rdata")

# Calculate seurat_cluster markers, with high stringency to return very
# specific markers, quickly.
seurat_cluster_markers <- FindAllMarkers(obj,
                                  assay = "RNA",
                                  logfc.threshold = 1.5,
                                  min.pct = 0.5,
                                  only.pos = T,
                                  return.thresh = 0.0001,
                                  densify = T,
                                  verbose = T)
write.csv(seurat_cluster_markers,
          file = "results/basic-annotation/seurat_cluster_markers.csv",
          row.names = F)
top5 <- seurat_cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# Read in canonical markers
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
  print(p1 + p2 + plot_layout(ncol = 2, widths = c(1, 2)))
  dev.off()
}

# Export a basic DimPlot of clusters
p <- DimPlot(obj, label = T, raster = F)
pdf(file = "results/basic-annotation/DimPlot.pdf",
    useDingbats = F)
print(p)
dev.off()

###############################################################################
# Annotation notes
###############################################################################

# Macrophages
# 0 - Mac-1 - Fcgr1+ Adgre1+ Cd68+ Lgals3+
# 1 - Mac-2 - H2-Ab1+ H2-Aa+ Cd74+ Fcgr1+ Adgre1+ Cd68+ Lgals3+
# 2 - Mac-3 - Fcgr1+ Adgre1+ Cd68+ Lgals3+
# 6 - Mac-4 - Fcgr1+ Adgre1+ Cd68+ Lgals3+
# 10 - Mac-5 - Fcgr1+ Adgre1+ Cd68+ Lgals3+

# DC-like cells
# 13 - DC-1 - Cd209a+, CD11c+ (Itgax), H2-Ab1+ H2-Aa+ Cd74+ Fcgr1+ Cd68+ Lgals3+
# 18 - DC-2 - H2-Ab1+ H2-Aa+ Cd74+ Cd68+ Lgals3+
#      Likely conventional DC, Itgax+ Naaa+ Irf8+ Xcr1+ Clec9a+ PMID: 29925006
# 21 - DC-3 - H2-Ab1+ H2-Aa+ Cd74+ Ccl5+ (DC)
#      Likely migratory DC, Ccr7+ Fscn1+ PMID: 29925006

# Endothelial cells
# 3 - EC-1 - Cdh5+ Pecam1+ Kdr+ Fabp4+
# 11 - EC-2 - Cdh5+ Pecam1+ Kdr+ Fabp4+
# 19 - EC-3 - Cdh5+ Pecam1+ Kdr+ Fabp4+

# Fibroblasts
# 5 - Fibro-1 - Pdgfra+ Tcf21+ Col1a1+ Postn++ Cthrc1++
# 12 - Fibro-2 - Pdgfra+ Tcf21+ Col1a1+ Gsn++ Postn+
# 14 - Fibro-3 - Pdgfra+ Tcf21+ Col1a1+ Postn+ Cthrc1+ Mki67+ Ccnb2+ Ccn2a+
#      Stmn1+, Wt1+ Dmkn+ Saa3+ Krt8+ Krt19+, Mixture of activating,
#      proliferating fibroblasts and epicardial cells

#Cycling cells
# 7 - Cycling - Mki67+ Ccnb2+ Ccna2+ Stmn1++ Fcgr1+ Adgre1+ Cd68+ Lgals3+,
#     Proliferating and mostly macrophages, but not entirely
#     macrophages, there are some endothelial cells and fibroblasts inside

# Mural cells (smooth muscle cells and pericytes)
# 20 - Mural - Myh11+ Rgs5+ Tagln+ Pdgfrb+ Cspg4+ Vtn+ Postn+ Fabp4+

# Granulocytes
# 4 - Gran-1 - S100a8+ S100a9+ Csf3r+
# 9 - Gran-2 - S100a8+ S100a9+ Csf3r+

# B cells
# 8 - B-cell - Ms4a1+ Cd79a+ H2-Ab1+ Cd74+  H2-Aa+

# T cells
# 15 - T-cell-1 - Cd3e+ Cd3d+ Ccl5+
# 16 - T-cell-2 - Cd3e+ Cd3d+ Lef1+

# NK cells
# 17 - NK-cell - Ncr1+ Klrk1+ Ccl5+

###############################################################################

# Rename Idents to annotations
obj <- RenameIdents(obj,
                    "0" = "Mac-1",
                    "1" = "Mac-2",
                    "2" = "Mac-3",
                    "3" = "EC-1",
                    "4" = "Gran-1",
                    "5" = "Fibro-1",
                    "6" = "Mac-4",
                    "7" = "Cycling",
                    "8" = "B-cell",
                    "9" = "Gran-2",
                    "10" = "Mac-5",
                    "11" = "EC-2",
                    "12" = "Fibro-2",
                    "13" = "DC-1",
                    "14" = "Fibro-3",
                    "15" = "T-cell-1",
                    "16" = "T-cell-2",
                    "17" = "NK-cell",
                    "18" = "DC-2",
                    "19" = "EC-3",
                    "20" = "Mural",
                    "21" = "DC-3"
                    )

# Store renamed idents as a new meta data column, set as Idents
obj@meta.data$basic_annotation <- Idents(obj)

# Refactor annotation levels
source("scripts/etc/dimplotlevels.R")
obj@meta.data$basic_annotation <- factor(obj@meta.data$basic_annotation,
                                         levels = dimplotlevels)
DimPlot(obj,
        group.by = "basic_annotation",
        label = T,
        repel = T)

# Set Idents as re-factored basic_annotation identities
Idents(obj) <- "basic_annotation"

# Save object with basic annotations
save(obj, file = "results/objects/obj_annotated.Rdata")

# Saved basic annotation, barcodes and UMAP embeddings etc. for consistency in
# external usage, de-comment to overwrite
barcodes <- rownames(obj@meta.data)
annotation <- obj@meta.data$basic_annotation
genotype <- obj@meta.data$genotype
timepoint <- obj@meta.data$timepoint
sample <- obj@meta.data$sample
umap_1 <- Embeddings(obj[["umap"]])[,1]
umap_2 <- Embeddings(obj[["umap"]])[,2]
basic_annotation <- data.frame(barcodes,
                               annotation,
                               genotype,
                               timepoint,
                               sample,
                               umap_1,
                               umap_2)
write.csv(basic_annotation, file = "data/basic_annotation.csv", row.names = F)

# Quality control summary of basic annotated object

# Counting cluster markers
marker_count <- seurat_cluster_markers %>%
  group_by(cluster) %>%
  summarise(nMarkers = n())
marker_count$cluster <- as.character(marker_count$cluster)

# Finding most common phase designation
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

# Summarizing QC metrics
qc_summary <- obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarize(mean(percent.mt),
            mean(nFeature_RNA),
            mean(nCount_RNA),
            calculate_mode(basic_annotation))
qc_summary <- qc_summary %>% left_join(marker_count,
                                       by = c("seurat_clusters" = "cluster"))

colnames(qc_summary) <- c("cluster",
                          "mean_percent_mt",
                          "mean_nFeatures",
                          "mean_nCounts",
                          "basic_annotation",
                          "nMarkers")
write.csv(qc_summary,
          file = "results/basic-annotation/final_qc_summary.csv",
          row.names = F)

# Plotting QC summary
p <- ggplot(qc_summary, aes(x = mean_percent_mt,
                            y = mean_nFeatures,
                            colour = nMarkers,
                            label = basic_annotation)) +
  scale_colour_viridis_c() +
  geom_point(size = 3) +
  geom_text_repel() +
  ggtitle("Final clustering QC")
pdf(file = "results/basic-annotation/final_qc_summary_scatter.pdf",
    width = 8,
    height = 6,
    useDingbats = F)
print(p)
dev.off()

# Plotting QC overview summary + nFeature + Clusters
p1 <- DimPlot(obj, label = T, raster = F) + NoLegend()
p2 <- FeaturePlot(obj,
                  features = "nFeature_RNA",
                  label = T,
                  raster = F) + NoLegend()
p3 <- ggplot(qc_summary, aes(x = mean_percent_mt,
                             y = mean_nFeatures,
                             colour = nMarkers,
                             label = basic_annotation)) +
  scale_colour_viridis_c() +
  geom_point(size = 3) +
  geom_text_repel() +
  ggtitle("Final clustering QC")
pdf(file = "results/basic-annotation/final_qc_overview.pdf",
    width = 13,
    height = 10,
    useDingbats = F)
print((p1 / p2) | p3 )
dev.off()
