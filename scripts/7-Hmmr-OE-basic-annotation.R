# Basic cluster annotation based on canonical cell-type markers

# Load libraries
library(Seurat) #v4.0.1
library(patchwork)
library(ggplot2)
library(dplyr)

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
# 1 - Mac-2 - Fcgr1+ Adgre1+ Cd68+ Lgals3+
# 2 - Mac-3 - H2-Ab1+ H2-Aa+ Cd74+ Fcgr1+ Adgre1+ Cd68+ Lgals3+
# 9 - Mac-4 - Cd68+ Lgals3+
# 11 - Mac-5 - Fcgr1+ Adgre1+ Cd68+ Lgals3+
# 13 - Mac-6 - Fcgr1+ Adgre1+ Cd68+ Lgals3+
# 14 - Mac-7 - H2-Ab1+ H2-Aa+ Cd74+ Fcgr1+ Adgre1+ Cd68+ Lgals3+
# 16 - Mac-8 - Fcgr1+ Adgre1+ Cd68+ Lgals3+
# 21 - Mac-9 - Fcgr1+ Adgre1+ Cd68+ Lgals3+

# DC-like cells
# 15 - DC-1 - Cd209a+, CD11c+ (Itgax), H2-Ab1+ H2-Aa+ Cd74+ Fcgr1+ Cd68+ Lgals3+
# 24 - DC-2 - H2-Ab1+ H2-Aa+ Cd74+ Cd68+ Lgals3+
#      Likely conventional DC, Itgax+ Naaa+ Irf8+ Xcr1+ Clec9a+ PMID: 29925006
# 30 - DC-3 - H2-Ab1+ H2-Aa+ Cd74+ Ccl5+ (DC)
#      Likely migratory DC, Ccr7+ Fscn1+ PMID: 29925006

# Endothelial cells
# 3 - EC-1 - Cdh5+ Pecam1+ Kdr+ Fabp4+
# 10 - EC-2 - Cdh5+ Pecam1+ Kdr+ Fabp4+
# 17 - EC-3 - Cdh5+ Pecam1+ Kdr+ Fabp4+
# 19 - EC-4 - Cdh5+ Pecam1+ Kdr+ Fabp4+
# 25 - EC-5 - Cdh5+ Pecam1+ Kdr+ Fabp4+
# 28 - EC-6 - Cdh5+ Pecam1+ Kdr+ Fabp4+
# 33 - EC-7 - Pecam1+ Fabp4+

# Fibroblasts
# 4 - Fibro-Myo - Pdgfra+ Tcf21+ Col1a1+ Postn++ Cthrc1++
# 12 - Fibro-Rest - Pdgfra+ Tcf21+ Col1a1+ Gsn++ Postn+
# 20 - Fibro-Act - Col1a1+ Postn+ Cthrc1+

# Epicardial cells
# 26 - Epi - Wt1+ Dmkn+ Saa3+ Krt8+ Krt19+

# Mural cells (smooth muscle cells and pericytes)
# 29 - Mural - Myh11+ Rgs5+ Tagln+ Pdgfrb+ Cspg4+ Vtn+ Postn+ Fabp4+

# Granulocytes
# 5 - Gran-1 - S100a8+ S100a9+ Csf3r+
# 6 - Gran-2 - S100a8+ S100a9+ Csf3r+
# 31 - Gran-3 - S100a8+ S100a9+ Csf3r+

# B cells
# 7 - B-cell - Ms4a1+ Cd79a+ H2-Ab1+ Cd74+  H2-Aa+

# T cells
# 18 - T-cell-1 - Cd3e+ Cd3d+ Ccl5+
# 22 - T-cell-2 - TCd3e+ Cd3d+ Lef1+

# NK cells
# 23 - NK-cell - Ncr1+ Klrk1+ Ccl5+

# Cardiomyocytes
# 27 - CM - Actc1+ Nppa+ Nppb+ Tnnc1+ Tnnt2+

# Cycling
# 8 - Cyc-1 - Mki67+ Ccnb2+ Ccna2+ Stmn1++ Fcgr1+ Adgre1+ Cd68+ Lgals3+,
#     Cluster 8 is proliferating and mostly macrophages, but not entirely
#     macrophages, there are some endothelial cells and fibroblasts inside
# 32 - Cyc-2 - Mki67+ Ccnb2+ Ccna2+ Stmn1++ Cd3e+ Cd3d+ Ccl5+,
#      Cluster 32 is mostly proliferating T-cells, but to be safe I will,
#      annotate these as the second cycling cluster
###############################################################################

# Rename Idents to annotations
obj <- RenameIdents(obj,
                    "0" = "Mac-1",
                    "1" = "Mac-2",
                    "2" = "Mac-3",
                    "3" = "EC-1",
                    "4" = "Fibro-Myo",
                    "5" = "Gran-1",
                    "6" = "Gran-2",
                    "7" = "B-cell",
                    "8" = "Cycling-1",
                    "9" = "Mac-4",
                    "10" = "EC-2",
                    "11" = "Mac-5",
                    "12" = "Fibro-Rest",
                    "13" = "Mac-6",
                    "14" = "Mac-7",
                    "15" = "DC-1",
                    "16" = "Mac-8",
                    "17" = "EC-3",
                    "18" = "T-cell-1",
                    "19" = "EC-4",
                    "20" = "Fibro-Act",
                    "21" = "Mac-9",
                    "22" = "T-cell-2",
                    "23" = "NK-cell",
                    "24" = "DC-2",
                    "25" = "EC-5",
                    "26" = "Epi",
                    "27" = "CM",
                    "28" = "EC-6",
                    "29" = "Mural",
                    "30" = "DC-3",
                    "31" = "Gran-3",
                    "32" = "Cycling-2",
                    "33" = "EC-7")

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
# barcodes <- rownames(obj@meta.data)
# annotation <- obj@meta.data$basic_annotation
# genotype <- obj@meta.data$genotype
# timepoint <- obj@meta.data$timepoint
# sample <- obj@meta.data$sample
# UMAP_1 <- Embeddings(obj[["umap"]])[,1]
# UMAP_2 <- Embeddings(obj[["umap"]])[,2]
# basic_annotation <- data.frame(barcodes,
#                                annotation,
#                                genotype,
#                                timepoint,
#                                sample,
#                                UMAP_1,
#                                UMAP_2)
# write.csv(basic_annotation, file = "data/basic_annotation.csv", row.names = F)
