# Differential gene expression analysis of all interstitial cells

# Load libraries
library(Seurat) #>=4.0.1
library(ggplot2)
library(dplyr)
source("scripts/etc/dimplotlevels.R")

# Load object
load("results/objects/obj_annotated.Rdata")

# Seurat based Wilcox method, no threshold, loop through each cluster
genes <- list()
for (i in levels(Idents(obj))) {
  results <- FindMarkers(obj,
                         subset.ident = i,
                         group.by = "genotype",
                         ident.1 = "OE",
                         base = 2,
                         logfc.threshold = 0,
                         densify = T)
  results$cluster <- i
  results$gene <- row.names(results)
  row.names(results) <- NULL
  results$regulation <- ifelse(results$avg_log2FC > 0, "Up", "Down")
  genes[[i]] <- results
}
dge_no_threshold <- do.call(rbind, genes)
rownames(dge_no_threshold) <- NULL
write.csv(dge_no_threshold,
          file = "results/differential-gene-expression/dge_no_threshold.csv",
          row.names = F)

# Filter out non-significant genes
dge <- dge_no_threshold[dge_no_threshold$p_val_adj < 0.01 &
                          abs(dge_no_threshold$avg_log2FC) > 0.25, ]

# Save significant deg
write.csv(dge,
          file = "results/differential-gene-expression/dge.csv",
          row.names = F)

# Count and plot deg per cluster
dge$regulation <- factor(dge$regulation, levels = c("Up", "Down"))
dge$cluster <- factor(dge$cluster, levels = rev(dimplotlevels))
p <- dge %>%
  group_by(cluster, regulation) %>%
  ggplot(aes(x = cluster)) +
  geom_bar(aes(fill = regulation)) +
  coord_flip() +
  scale_fill_manual(name = "Regulation in OE",
    values = c("#D75438", "#4878CD")) +
  theme_bw() +
  ggtitle("Number of differentially expressed genes")
pdf(file = "results/differential-gene-expression/dge_count.pdf",
    useDingbats = F)
print(p)
dev.off()
