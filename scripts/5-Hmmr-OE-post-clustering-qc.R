# Load libraries
library(Seurat) #v4.0.1
library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(ggrepel)

# Load object
load("results/objects/obj.Rdata")

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
p1 <- FeaturePlot(obj, features="diss_genes1", label = T, raster = F) + 
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
p1 <- FeaturePlot(obj, features="percent.mt", label = T, raster = F)
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
g2m_genes_no_Hmmr <- g2m_genes[-42]
remove(cc_genes) 
obj <- CellCycleScoring(obj,
                        s.features = s_genes,
                        g2m.features = g2m_genes_no_Hmmr)
obj@meta.data$Phase <- factor(obj@meta.data$Phase,
                              levels = c("G1", "S", "G2M"))
phase_membership <- (prop.table(table(obj$Phase, Idents(obj)),
                                margin = 2)*100)
phase_membership <- as.data.frame(phase_membership)
colnames(phase_membership) <- c("Phase", "Cluster", "Percent")
phase_membership$Percent.round <- round(phase_membership$Percent)
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
  geom_text(aes(label = paste0(Percent.round, "%", sep = "")),
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
write.csv(top5,
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

###############################################################################
# Manual investigation of clusters that may be low-quality cells.
###############################################################################

# Investigating clusters with <1000 mean_nfeatures, cluster numbering might be 
# different based on your machine, for me, these are specifically clusters:
# 10, 23, 29 and 33

# Cluster 10: Likely macrophage sub-population, percent_mt is fine, G1, 
# low number marker genes (116) but they seem to have plausible macrophage
# function.

# Cluster 23: Likely macrophage sub-population, percent_mt is high, G2M, 
# low number of markers genes (86), but function is consistent with macrophages,
# could be macs that are starting or ending proliferation.

# Cluster 29: Cardiomyocytes, expressing canonical cardiomyocyte marker genes,
# Actc1, Tnnt2, but also express Ankrd1 a known border zone gene. Perhaps these
# are border zone cardiomyocytes small enough to make into through FACS 
# and into gems.

# Cluster 33: Likely macrophage sub-population, percent_mt is ok, G1, lowest
# number of marker genes. Top markers are lncRNA (Gm42418, Gm26917), 
# interestingly other groups have found and excluded clusters with these 
# top 2 markers genes (PMID: 33205009). Also has stress markers Jun, Fosb. 
# Likely low-quality macrophages.

# Investigating clusters with >5% mean_percent_mt, cluster numbering might be 
# different based on your machine, for me, these are specifically clusters:
# 8, 9, 23, 35

# Cluster 8: B-cells, with typical low nFeature. Express canonical
# B-cell markers Igkc, Cd79a, Ms4a1.

# Cluster 9: Endothelial sub-population, G1, 146 marker genes, top markers
# have plausible EC function (Fabp4, Cd36)

# Cluster 23: See above.

# Cluster 35: Likley endothelial cell sub-population, G1, 415 marker genes.
# Top markers have plausible EC function, including Plvap, 
# an EC-specific membrane protein.

###############################################################################

# Sub-setting object to exclude low-quality clusters
obj <- subset(x = obj, idents = "33", invert = TRUE)

# Remove extra meta data columns
obj@meta.data <- select(obj@meta.data, -starts_with("pANN"), -c(nCount_SCT,
                                              nFeature_SCT,
                                              SCT_snn_res.0.8,
                                              seurat_clusters,
                                              doublet_classification),)
# Keep only necessary counts and data slots 
obj <- DietSeurat(obj,
                  counts = T,
                  data = T,
                  scale.data = F,
                  features = NULL,
                  assays = NULL,
                  dimreducs = NULL,
                  graphs = NULL)

# Save object
save(obj, file = "results/objects/obj.Rdata")
