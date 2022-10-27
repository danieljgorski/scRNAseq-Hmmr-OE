# Differential progression analysis on fibroblasts, workflow adopted from 
# https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html
# and https://hectorrdb.github.io/condimentsPaper/articles/TGFB.html

# Load libraries
library(Seurat) #v4.0.1
library(ggplot2)
library(patchwork)
library(condiments)
library(tradeSeq)
library(slingshot)
library(dplyr)
library(SingleCellExperiment)
library(RColorBrewer)
library(cowplot)
library(scales)
library(pheatmap)
library(scater)
source("scripts/etc/geno_colors.R")

# Load the fibroblast subset object
load("results/objects/fibro.Rdata")

# Remove epicardial and cycling cells
fibro <- subset(x = fibro,
                subset = Mki67 > 0.1 | Ccnb2 > 0 | Wt1 > 0 | Clu > 0,
                invert = T)
DimPlot(fibro, label = T)

# Convert to singleCellExperiment
sce <- as.SingleCellExperiment(fibro, assay = "RNA")

# Plot the genotypes
df <- bind_cols(
  as.data.frame(reducedDims(sce)$UMAP),
  as.data.frame(colData(sce)[, -3])
) %>%
  sample_frac(1)
p1 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = genotype)) +
  geom_point(size = .7) +
  scale_color_manual(values=c("#999999", "#008000")) +
  labs(col = "Genotype", x="UMAP-1", y="UMAP-2") +
  theme_classic() +
  theme(axis.text = element_blank())
p1
pdf(file = "results/fibro-differential-progression/genotypes.pdf",
    height = 4,
    width = 6,
    useDingbats = F)
print(p1)
dev.off()

# Calculate the imbalance score and visualize
scores <- condiments::imbalance_score(
  Object = df %>% select(UMAP_1, UMAP_2) %>% as.matrix(), 
  conditions = df$genotype,
  k = 20, smooth = 40)
df$scores <- scores$scaled_scores
p3 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = scores)) +
  geom_point(size = .7) +
  scale_color_viridis_c(option = "C") +
  labs(col = "Imbalance\nscore", x="UMAP-1", y="UMAP-2") +
  theme_classic() +
  theme(axis.text = element_blank())
p3
pdf(file = "results/fibro-differential-progression/imbalance.pdf",
    height = 4,
    width = 6,
    useDingbats = F)
print(p3)
dev.off()

# Trajectory inference with slingshot, starting in Fibro-2, which are
# are resting fibroblasts (Gsn+), ending in Fibro-1, which are activated
# fibroblasts, (Postn+, Cthrc1+)
sce <- slingshot(sce, reducedDim = 'UMAP',
                 clusterLabels = colData(sce)$basic_annotation,
                 start.clus = 'Fibro-2',
                 end.clus = "Fibro-1")

# Test whether trajectory should be fitted independently for
# different conditions or not
set.seed(821)
topologyTest(SlingshotDataSet(sce),
             sce$genotype,
             rep = 100,
             methods = "KS_mean",
             threshs = .01)

# Plot a joint trajectory, topology test p-value = 1, not significant
df <- bind_cols(
  as.data.frame(reducedDims(sce)$UMAP),
  as.data.frame(colData(sce)[,-16]) # cannot include the pseudotime ordering column here
) %>%
  sample_frac(1)
curve <- slingCurves(sce)[[1]]
p4 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = slingPseudotime_1)) +
  geom_point(size = .7) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime", x="UMAP-1", y="UMAP-2") +
  geom_path(data = curve$s[curve$ord, ] %>% as.data.frame(),
            col = "black", size = 1.25) +
  theme_classic() +
  theme(axis.text = element_blank())
p4
pdf(file = "results/fibro-differential-progression/trajectory.pdf",
    height = 4,
    width = 6,
    useDingbats = F)
print(p4)
dev.off()

# Plot the pseudotime distributions of each genotype
p5 <- ggplot(df, aes(x = slingPseudotime_1,
                     fill = genotype,
                     color = genotype)) +
  geom_density(alpha = .5) +
  scale_fill_manual(values=c("#999999", "#008000")) +
  scale_color_manual(values=c("#999999", "#008000")) +
  labs(x = "Pseudotime", fill = "Genotype", y="Density", color = "Genotype") +
  theme_classic() +
  theme(legend.position = "top",
        legend.justification = "left")
p5
pdf(file = "results/fibro-differential-progression/differential-progression.pdf",
    height = 4,
    width = 6,
    useDingbats = F)
print(p5)
dev.off()

# Kolmogorov-Smirnov Test for differential progression
progressionTest(SlingshotDataSet(sce), conditions = sce$genotype)
# KS-test statistic = 0.0877, p-value = 6.66e-16, significant

# Differential gene expression
set.seed(3)
icMat <- evaluateK(counts = sce,
                   conditions = factor(colData(sce)$genotype),
                   nGenes = 300,
                   k = 3:7)

# Fit GAM, 4 knots had best fit
set.seed(3)
sce <- fitGAM(counts = sce,
              nknots = 4,
              conditions = factor(colData(sce)$genotype))
mean(rowData(sce)$tradeSeq$converged)

# Differential expression across conditions through pseudotime
condRes <- conditionTest(sce, l2fc = log2(1.25))
condRes$padj <- p.adjust(condRes$pvalue, "fdr")
condRes <- condRes %>% filter(padj <= 0.05) %>% arrange(desc(waldStat))
high_exp <- assay(sce, "counts") %>% log1p() %>% rowMeans() # calculate log1p mean expression of each gene 
high_exp <- high_exp[high_exp > 0.4] # genes with high expression
high_exp <- names(high_exp)
condRes <- condRes %>% filter(rownames(condRes) %in% high_exp) # Filter out differentially expressed genes with low expression
write.csv(condRes, file = "results/fibro-differential-progression/deg_ps.csv")

# Heatmaps of genes DE between conditions, ordered according to a hierarchical 
# clustering on the WT condition
yhatSmooth <- predictSmooth(sce,
                            gene = rownames(condRes),
                            nPoints = 100,
                            tidy = FALSE) %>%
  log1p()
yhatSmoothScaled <- t(apply(yhatSmooth,1, scales::rescale))
heatSmooth_wt <- pheatmap(yhatSmoothScaled[, 1:100],
                          cluster_cols = FALSE,
                          border_color = NA,
                          show_rownames = FALSE,
                          show_colnames = FALSE,
                          main = "WT",
                          legend = FALSE,
                          silent = TRUE)
matchingHeatmap_oe <- pheatmap(yhatSmoothScaled[heatSmooth_wt$tree_row$order, 101:200],
                               cluster_cols = FALSE,
                               border_color = NA,
                               cluster_rows = FALSE,
                               show_rownames = TRUE,
                               show_colnames = FALSE,
                               main = "Hmmr-OE",
                               legend = FALSE,
                               silent = TRUE,
                               fontsize_row = 6)
p9 <- plot_grid(heatSmooth_wt[[4]], matchingHeatmap_oe[[4]], ncol = 2)
p9
pdf(file = "results/fibro-differential-progression/deg_ps_heatmap.pdf",
    height = 6,
    width = 6,
    useDingbats = F)
print(p9)
dev.off()

# Visualize genes
for (i in rownames(condRes)) {
  p <- plotSmoothers(sce,
                     assays(sce)$counts,
                     gene = i,
                     lwd = 2.5,
                     border = TRUE,
                     size = .3,
                     curvesCols = geno_colors) +
    scale_color_manual(values = geno_colors,
                       labels = c("WT", "Hmmr-OE")) +
    ggtitle(i) +
    theme_classic() +
    theme(plot.title = element_text(face = "italic", size = 24),
          legend.position = "none",
          legend.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 3.5)))
  pdf(file = paste0("results/fibro-differential-progression/deg_ps_", i, ".pdf"),
      height = 4,
      width = 6,
      useDingbats = F)
  print(p)
  dev.off()
}

# Save fibroblast slingshot object
save(sce, file = "results/objects/fibro-slingshot.Rdata")