# Plotting common differentially downregulated genes in macrophage clusters,
# and highlighting genes involved in cell migration and leukocyte chemotaxis

# Load libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(forcats)
source("scripts/etc/maclevels.R")

# Read in dge and cell_migration GOBP data
dge_full <- read.csv("results/mac-differential-gene-expression/mac_dge_no_threshold.csv")
cell_migration <- read.csv("data/cell_migration_GO_0016477.csv")
cell_migration <- unique(cell_migration$Symbol)
leukocyte_chemotaxis <- read.csv("data/leukocyte_chemotaxis_GO_0030595.csv")
leukocyte_chemotaxis <- unique(leukocyte_chemotaxis$Symbol)

# Select significantly down regulated genes in more than one cluster
dge <- dge_full %>%
  filter(abs(avg_log2FC) > 0.25) %>%
  filter(p_val_adj < 0.05) %>%
  filter(regulation == "Down") %>%
  group_by(gene) %>%
  filter(n() > 1)

# Summarize mean pvalue, avg_logFC, number of clusters, arrange by number of
# clusters with regulation, then logFC, then p_val_adj
df <- dge %>%
  group_by(gene) %>%
  summarize(mean_p_val_adj = mean(p_val_adj),
            mean_avg_log2FC = mean(avg_log2FC),
            count = n_distinct(cluster)) %>%
  arrange(desc(count), mean_avg_log2FC, mean_p_val_adj)

# Factor genes for consistent plot ordering
df$gene <- factor(df$gene, levels = rev(c(df$gene)))

# Find if gene is listed in cell migration GO term, if so assign color and
# GO term
df$GO <- "maybe"
for (row in 1:nrow(df)) {
  if (df[row, ]$gene %in% cell_migration) {
    df[row, ]$GO <- "migration"
    if (df[row, ]$gene %in% leukocyte_chemotaxis) {
      df[row, ]$GO <- "chemotaxis"
    }
  } else{
    df[row, ]$GO <- "other"
  }
}


# Plot
pdf(file = "results/mac-common-dge/mac-common-down-dge.pdf",
    useDingbats = F)
ggplot(df, aes(x = count, y = gene)) +
  geom_segment(
    aes(x = 0, xend = count, y = gene, yend = gene),
    color = ifelse(df$GO %in% c("migration", "chemotaxis"), "#0458AB", "gray70"),
    size = ifelse(df$GO %in% c("migration", "chemotaxis"), 0.8, 0.5)) +
  geom_point(
    color = ifelse(df$GO %in% c("migration", "chemotaxis"), "#0458AB", "gray70"),
    size=ifelse(df$GO %in% c("migration", "chemotaxis"), 1.7, 1.2)) +
  geom_text(
    aes(label = gene,
        hjust = 0,
        fontface = "italic"),
    position = position_nudge(x = 0.2),
    color = "#0458AB",
    size = ifelse(df$GO == "chemotaxis", 3.25, 0)) +

  ggtitle(expression(paste("Reduced chemotaxis in ",
                           italic("Hmmr"),
                           " OE macrophages"))) +
  ylab("Downregulated in > 1 macrophage cluster") +
  xlab("Number of macrophage clusters with downregulation") +

  scale_x_continuous(limits = c(0, 13),
                     breaks = c(0, 2, 4, 6, 8, 10, 12),
                     expand = c(0, 0.01)) +
  coord_cartesian(clip = "off") +

  theme_minimal() +
  theme(plot.title = element_text(color = "#5c5c5c", size = 14),
        plot.margin = unit(c(1, 3, 1, 1), "cm"),
        panel.grid = element_blank(),
        axis.text.y = element_text(color = "gray70", size = 8),
        axis.title.y = element_text(hjust = 0.5, color = "gray70", size = 10),
        axis.title.x = element_text(hjust = 0, color = "gray70", size = 10),
        axis.text.x = element_text(color = "gray70", size = 8)) +

  annotate("text", x = 7.5, y = 15, label = "cell migration GO:0016477",
           color = "#0458AB", size = 3.5, hjust = 0) +
  annotate("segment",
           x = 7.95,
           xend = 7.95,
           y = 14,
           yend = 13,
           colour = "#0458AB",
           size = 0.25) +
  annotate("segment",
           x = 7.95,
           xend = 8.25,
           y = 13,
           yend = 13,
           colour = "#0458AB",
           size = 0.25) +
  annotate("text", x = 8.35, y = 13, label = "leukocyte chemotaxis GO:0030595",
           color = "#0458AB", size = 3.5, fontface = "italic", hjust = 0) +
  annotate("text", x = 7.5, y = 11.25, label = "other",
           color = "gray70", size = 3.5, hjust = 0)
dev.off()

# Write df
write.csv(df,
          file = "results/mac-common-dge/mac-common-down-dge.csv",
          row.names = F)

# Find number of downregulated leukocyte chemotaxis in each cluster
df <- dge %>% group_by(cluster) %>%
  tally(gene %in% leukocyte_chemotaxis) %>%
  arrange(desc(n))
df$cluster <- factor(df$cluster, levels = rev(df$cluster))

# Plot
pdf(file = "results/mac-common-dge/leukocyte_chemotaxis_profile.pdf",
    useDingbats = F)
ggplot(df, aes(x = n, y = cluster)) +
  geom_bar(stat = "identity", fill = "#0458AB") +
  scale_x_continuous(limits = c(0, 10),
                    breaks = c(0, 2, 4, 6, 8, 10),
                    expand = c(0.01, 0.01)) +
  ggtitle("Number of downregulated chemotaxis genes per cluster") +
  xlab("Number of downregulated genes involved in leukocyte chemotaxis") +
  theme_void() +
  theme(plot.title = element_text(color = "#5c5c5c", size = 14),
        plot.margin = unit(c(1, 2, 1, 1), "cm"),
        axis.text.y = element_text(colour = "gray70", size = 10, hjust = 1),
        axis.text.x = element_text(colour = "gray70", size = 10),
        axis.title.x = element_text(colour = "gray70", size = 11, hjust = 0, vjust = -1),
        legend.position = c(0.95, 0.15),
        legend.key.height = unit(.3, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.title = element_text(colour = "gray70"),
        legend.title.align = -2,
        legend.text = element_text(colour = "gray70"))
dev.off()
