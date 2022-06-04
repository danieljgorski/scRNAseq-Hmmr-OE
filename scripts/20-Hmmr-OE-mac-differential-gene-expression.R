# Differential gene expression analysis of macrophage subset

# Load libraries
library(Seurat) #>=4.0.1
library(ggplot2)
library(dplyr)
source("scripts/etc/maclevels.R")

# Load object
load("results/objects/mac_annotated.Rdata")

# Seurat based Wilcox method, no threshold, loop through each cluster
genes <- list()
for (i in levels(Idents(mac))) {
  results <- FindMarkers(mac,
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
mac_dge_no_threshold <- do.call(rbind, genes)
rownames(mac_dge_no_threshold) <- NULL
write.csv(mac_dge_no_threshold,
          file = "results/mac-differential-gene-expression/mac_dge_no_threshold.csv",
          row.names = F)

# Filter out non-significant genes
mac_dge <- mac_dge_no_threshold[mac_dge_no_threshold$p_val_adj < 0.01 &
                          abs(mac_dge_no_threshold$avg_log2FC) > 0.25, ]

# Save significant deg
write.csv(mac_dge,
          file = "results/mac-differential-gene-expression/mac_dge.csv",
          row.names = F)

# Count and plot deg per cluster
mac_dge$regulation <- factor(mac_dge$regulation, levels = c("Up", "Down"))
mac_dge$cluster <- factor(mac_dge$cluster, levels = rev(maclevels))
p <- mac_dge %>%
  group_by(cluster, regulation) %>%
  ggplot(aes(x = cluster)) +
  geom_bar(aes(fill = regulation)) +
  coord_flip() +
  scale_fill_manual(name = "Regulation in OE",
                    values = c("#D75438", "#4878CD")) +
  theme_bw() +
  ggtitle("Number of differentially expressed genes")
pdf(file = "results/mac-differential-gene-expression/mac_dge_count.pdf",
    useDingbats = F)
print(p)
dev.off()
