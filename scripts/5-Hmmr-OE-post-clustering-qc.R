# Quality control after initial clustering

# Load libraries
library(Seurat) #v4.0.1
library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(ggrepel)

# Load object
load("results/objects/obj_integrated.Rdata")

# Switch default assay to RNA
DefaultAssay(obj) <- "RNA"

# Dissociation related genes
diss_genes <- unique(c("Atf3", "Btg2", "Cebpb", "Cebpb",
                       "Cxcl3", "Cxcl2", "Cxcl1",
                       "Dnaja1", "Dnajb1", "Dusp1",
                       "Egr1", "Fos", "Fosb", "Hsp90aa1",
                       "Hsp90ab1", "Hspa1a", "Hspa1b",
                       "Hspa1a", "Hspa1b", "Hspa8",
                       "Hspb1", "Hspe1", "Hsph1", "Id3",
                       "Ier2", "Jun", "Junb", "Jund",
                       "Mt1", "Nfkbia", "Nr4a1", "Ppp1r15a",
                       "Socs3", "Zfp36"))
obj <- AddModuleScore(obj,
                      features = list(diss_genes),
                      ctrl = 50,
                      name = "diss_genes")
p1 <- FeaturePlot(obj, features = "diss_genes1", label = T, raster = F) +
  ggtitle("Dissociation gene score")
p2 <- VlnPlot(obj, features = "diss_genes1", pt.size = 0, sort = T) +
  ggtitle("Dissociation gene score") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  NoLegend()
pdf(file = "results/post-clustering-qc/dissocation_genes.pdf",
    width = 14,
    height = 6,
    useDingbats = F)
p1 + p2 + plot_layout(ncol = 2)
dev.off()

# Batch effect exploration
p <- DimPlot(obj, group.by = "sample", raster = F)
pdf(file = "results/post-clustering-qc/DimPlot_sample.pdf")
print(p)
dev.off()

p <- DimPlot(obj, group.by = "date", raster = F)
pdf(file = "results/post-clustering-qc/DimPlot_date.pdf")
print(p)
dev.off()

p <- DimPlot(obj, group.by = "timepoint", raster = F)
pdf(file = "results/post-clustering-qc/DimPlot_timepoint.pdf")
print(p)
dev.off()

# nFeature_RNA
p1 <- FeaturePlot(obj, features = "nFeature_RNA", label = T, raster = F)
p2 <- VlnPlot(obj, features = "nFeature_RNA", pt.size = 0, sort = T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  NoLegend()
pdf(file = "results/post-clustering-qc/nFeature_RNA.pdf",
    width = 14,
    height = 6,
    useDingbats = F)
p1 + p2 + plot_layout(ncol = 2)
dev.off()

# nCount_RNA
p1 <- FeaturePlot(obj, features = "nCount_RNA", label = T, raster = F)
p2 <- VlnPlot(obj, features = "nCount_RNA", pt.size = 0, sort = T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  NoLegend()
pdf(file = "results/post-clustering-qc/nCount_RNA.pdf",
    width = 14,
    height = 6,
    useDingbats = F)
p1 + p2 + plot_layout(ncol = 2)
dev.off()

# percent.mt
p1 <- FeaturePlot(obj, features = "percent.mt", label = T, raster = F)
p2 <- VlnPlot(obj, features = "percent.mt", pt.size = 0, sort = T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  NoLegend()
pdf(file = "results/post-clustering-qc/percent_mt.pdf",
    width = 14,
    height = 6,
    useDingbats = F)
p1 + p2 + plot_layout(ncol = 2)
dev.off()

# Cell-cycle phase
cc_genes <- read_csv("data/seurat_cell_cycle.csv")
s_genes <- as.character(cc_genes$mouse.s.genes)
g2m_genes <- as.character(cc_genes$mouse.g2m.genes)
g2m_genes[42]
g2m_genes_no_hmmr <- g2m_genes[-42]
remove(cc_genes)
obj <- CellCycleScoring(obj,
                        s.features = s_genes,
                        g2m.features = g2m_genes_no_hmmr)
obj@meta.data$Phase <- factor(obj@meta.data$Phase,
                              levels = c("G1", "S", "G2M"))
phase_membership <- (prop.table(table(obj$Phase, Idents(obj)),
                                margin = 2) * 100)
phase_membership <- as.data.frame(phase_membership)
colnames(phase_membership) <- c("Phase", "Cluster", "Percent")
phase_membership$percent_round <- round(phase_membership$Percent)
write.csv(phase_membership,
          file = "results/post-clustering-qc/phase_membership.csv",
          row.names = F)
phase_membership$Phase <- factor(phase_membership$Phase,
                                 levels = c("G1", "S", "G2M"))
p1 <- DimPlot(obj,
              group.by = "Phase",
              cols = c("#e5e5e5", "#3a86ff", "#ffaa00"),
              raster = F)
p2 <- ggplot(phase_membership,
             aes(fill = Phase,
                 y = Percent,
                 x = Cluster)) +
  geom_bar(position = "stack", stat = "identity") +
  ylab("Percent of total") +
  xlab("Identity") +
  geom_text(aes(label = paste0(percent_round, "%", sep = "")),
            position = position_stack(vjust = 0.5), size = 2.5) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "Black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_manual(values = c("#e5e5e5", "#3a86ff", "#ffaa00"))
pdf(file = "results/post-clustering-qc/cell_cycle_phase.pdf",
    width = 17,
    height = 6,
    useDingbats = F)
p1 + p2 + plot_layout(guides = "collect")
dev.off()

# Find markers of initial clustering
cluster_markers <- FindAllMarkers(obj,
                                  assay = "RNA",
                                  logfc.threshold = 0.5,
                                  min.pct = 0.5,
                                  only.pos = T,
                                  return.thresh = 0.001,
                                  densify = T)
write.csv(cluster_markers,
          file = "results/post-clustering-qc/cluster_markers.csv",
          row.names = F)
top5 <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
write.csv(top5,
          file = "results/post-clustering-qc/cluster_markers_top_5.csv",
          row.names = F)
top10 <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
write.csv(top10,
          file = "results/post-clustering-qc/cluster_markers_top_10.csv",
          row.names = F)

# Counting cluster markers
marker_count <- cluster_markers %>%
  group_by(cluster) %>%
  summarise(nMarkers = n())
marker_count$cluster <- as.character(marker_count$cluster)

# Summarizing quality control
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}
qc_summary <- obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarize(mean(percent.mt),
            mean(diss_genes1),
            mean(nFeature_RNA),
            mean(nCount_RNA),
            calculate_mode(Phase))
qc_summary <- qc_summary %>% left_join(marker_count,
                                      by = c("seurat_clusters" = "cluster"))

colnames(qc_summary) <- c("cluster",
                          "mean_percent_mt",
                          "mean_dissociation_score",
                          "mean_nFeatures",
                          "mean_nCounts",
                          "cell_cycle",
                          "nMarkers")
write.csv(qc_summary,
          file = "results/post-clustering-qc/qc_summary.csv",
          row.names = F)

# Plotting QC summary
p <- ggplot(qc_summary, aes(x = mean_percent_mt,
                       y = mean_nFeatures,
                       colour = nMarkers,
                       label = cluster)) +
  scale_colour_viridis_c() +
  geom_point(aes(shape = cell_cycle), size = 3) +
  geom_text_repel() +
  ggtitle("Initial clustering QC")
pdf(file = "results/post-clustering-qc/qc_summary_scatter.pdf",
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
                            label = cluster)) +
  scale_colour_viridis_c() +
  geom_point(aes(shape = cell_cycle), size = 3) +
  geom_text_repel() +
  ggtitle("Initial clustering QC")
pdf(file = "results/post-clustering-qc/qc_overview.pdf",
    width = 13,
    height = 10,
    useDingbats = F)
print((p1 / p2) | p3 )
dev.off()

# Simple annotation, to identify remaining possible heterotypic multiplets

# Read in markers
markers <- read.csv(file = "data/canonical_markers.csv")
markers <- markers$All

# Loop through markers and generate Feature and VlnPlots of expression
for (i in markers) {
  p1 <- FeaturePlot(obj, features = i, label = T, raster = F)
  p2 <- VlnPlot(obj, features = i, pt.size = 0, sort = T) + NoLegend() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  pdf(file = paste0("results/post-clustering-qc/Feat_VlnPlot_", i, ".pdf"),
      width = 16,
      height = 5,
      useDingbats = F)
  print(p1 + p2 + plot_layout(ncol = 2, widths = c(1,2)))
  dev.off()
}

# Typical cell-types with top canonical markers

# Fibroblast markers, Col1a1, Tcf21, Postn
# Macrophage markers, Fcgr1, Adgre1,
# EC markers, Cdh5, Kdr, Pecam1
# Granulocytes, S100a8, S100a9,
# Proliferating, Mki67, Ccnb2
# DC_like Cd209a, H2-Ab1, Cd74
# Cardiomyocyte, Nppa, Actc1
# Bcell Ms4a1, Cd79a
# Tcell Cd3e, Cd3d
# NKcell Ncr1, Klrk1, Ccl5
# Mural Myh11, Vtn
# Epicardial Wt1, Dmkn, Krt19, Krt8
# Glial Plp1, Kcna1

# Notes on canonical marker expression across clusters

#0- Fcgr1, Adgre1
#1- Fcgr1, Adgre1
#2- Fcgr1, Adgre1, H2-Ab1, Cd74
#3- Cdh5, Kdr, Pecam1
#4- Ms4a1, Cd79a, H2-Ab1, Cd74
#5- S100a8, S100a9
#6- Fcgr1 low, Adgre1 low
#7- Fcgr1, Adgre1, Mik67, Ccnb2
#8- S100a8, S100a9
#9- Fcgr1, Adgre1, H2-Ab1 low, Cd74
#10- Fcgr1, Adgre1
#11- Cdh5, Kdr, Pecam1
#12- Cdh5 biphasic, Kdr biphasic, Pecam1 biphasic
#13- Fcgr1, Adgre1, H2-Ab1, Cd74
#14- Col1a1, Tcf21, Postn low
#15- Col1a1, Tcf21, Postn
#16- Fcgr1 low, Cd209a, H2-Ab1, Cd74, Klrk1 low
#17- Fcgr1, Adgre1 low, Cd74, Dmkn
#18- Col1a1, Tcf21, Postn
#19- Cd3e, Cd3d, Klrk1 low, Ccl5 low
#20- Cdh5, Kdr, Pecam1
#21- Col1a1 low, Tcf21 low, Postn low, Plp1 low, Kcna1 low
#22- Col1a1, Tcf21 low, Postn, Mik67 low, Ccnb2 low, Wt1, Dmkn
#23- Col1a1 low, Postn low, Fcgr1 low, Adgre1 low, S100a8 low, S100a9 low,
#    H2-Ab1, Cd74
#24- Fcgr1 low, Adgre1 low, H2-Ab1 low
#25- Fcgr1, Adgre1
#26- Cd3e, Cd3d
#27- S100a8, S100a9
#28- Ncr1, Klrk1, Ccl5
#29- H2-Ab1, Cd74, Klrk1
#30- Nppa, Nppb, Actc1
#31- Myh11, Vtn, Col1a1 low, Postn low
#32- Cdh5, Kdr, Pecam1, Wt1
#33- Fcgr1, Adgre1 low, Mik67, Ccnb2, H2-Ab1, Cd74
#34- Mik67 low, Ccnb2 low
#35- Cdh5, Kdr, Pecam1, Klrk1 low
#36- H2-Ab1 biphasic, Cd74 biphasic, Ccl5
#37- Cdh5 low, Pecam1 low, Cd209a, H2-Ab1, Cd74 biphasic
#38- Cdh5, Kdr, Pecam1, H2-Ab1, Cd74, Ms4a1, Cd79a
#39- Pecam1 low, S100a8 and S100a9 low, H2-Ab1, Cd74, Ms4a1 low, Cd79a biphasic
#40- Pecam1 low, H2-Ab1, Cd74, Ms4a1, Cd79a, Cd3e, Cd3d, Ccl5 low

# At risk clusters with multiple canonical cell type markers include:
# 23, 37, 38, 39, 40

###############################################################################
# Manual investigation of clusters that may be low-quality cells.
###############################################################################

# Investigating clusters with ~<1000 mean_nfeatures, cluster numbering might be
# different based on your machine, for me, these are specifically clusters:
# 5, 6, 12, 23, 24, 27, 30

# 5: Typical low feature granulocyte cluster with strong canonical
# marker gene expression (S100a8, S100a9). Low percent.mt

# 6: Seems to be macrophage cluster, but has low expression of Fcgr1, Adgre1
# and Cd68. But very high expression of Lgals3 (Mac-2). Marker genes include
# Fabp5, Ftl1, Fth1, Prdx1, all which plausible macrophage function. Low
# percent.mt

# 12: Contains cells with high and low endothelial marker gene expression, has
# relatively few marker genes (124), high percent.mt. However high expression
# of genes with plausible EC function, e.g. Fabp4 (known EC marker gene) and
# Cd36 (fatty acid transport).

# 23: Very high percent.mt, marker genes have a mix of cell types. Likely
# multiplets of fibroblasts, macrophages, granulocytes and DC-like cells.
# Will remove from the analysis.

# 24: Very few marker genes (58), likely macrophage but has low expression of
# Fcgr1 and Adgre1. High percent.mt. Marker genes include long-non-coding RNA
# (Gm42418 and Gm26917), similar to PMID: 33205009. tRNA synthetase Lars2,
# and transcription factor AY036118. Likely low-quality macrophages, will
# removed from analysis.

# 27: Typical low feature granulocyte cluster with strong canonical
# marker gene expression (S100a8, S100a9). Low percent.mt

# 30: Cardiomyocytes, strong expression of canonical markers Actc1, Nppa, Nppb.

# Investigating clusters with >5% mean_percent_mt, cluster numbering might be
# different based on your machine, for me, these are specifically clusters:
# 4, 12, 23, 24, 25

# 4: B-cells, strong expression of canonical markers Ms4a1, Cd79a

# 12: See above.

# 23: See above.

# 24: See above

# 25: Borderline high percent.mt, however has high nFeature, many markers genes
# (244) and strong expression of macrophage markers Fcgr1, Adgre1. Likely a
# macrophage sub-population.

# Investigating very small clusters, cluster numbering might be
# different based on your machine, for me, these are specifically clusters:
# 37, 38, 39, 40

# Unfortunately clusters 37, 38, 39, 40 all express a combination of endothelial
# cell, innate immune cell, and B-cell markers. Because its possible that these
# populations are multiplets of endothelial cells with bound extravasating
# leukocytes, they will be removed.

###############################################################################

# Exclude low-quality clusters determined during post-clustering-qc
obj <- subset(x = obj,
              idents = c("23",
                         "24",
                         "37",
                         "38",
                         "39",
                         "40"),
              invert = TRUE)

# Remove extra meta data columns created during integrated clustering
obj@meta.data <- select(obj@meta.data,
                        -c(nCount_SCT,
                           nFeature_SCT,
                           integrated_snn_res.0.8,
                           seurat_clusters))

# Keep only necessary counts and data slots
obj <- DietSeurat(obj,
                  counts = T,
                  data = T,
                  scale.data = F,
                  features = NULL,
                  assays = "RNA",
                  dimreducs = NULL,
                  graphs = NULL)

# Save clean object
save(obj, file = "results/objects/obj_clean.Rdata")

# Clear memory
rm(list = ls())
gc()
