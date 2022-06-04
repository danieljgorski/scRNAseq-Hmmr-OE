# Differential abundance analysis of macrophage subset

# Load libraries
library(Seurat) #>=4.0.1
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(scProportionTest)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(edgeR)
source("scripts/etc/colors.R")
source("scripts/etc/maclevels.R")

# Load object
load("results/objects/mac_annotated.Rdata")

# Proportion calculation----
mac_prop_by_sample <- (prop.table(table(mac@meta.data$mac_annotation,
                                          mac@meta.data$sample), margin = 2))
mac_prop_by_sample <- as.data.frame(mac_prop_by_sample)
colnames(mac_prop_by_sample) <- c("cluster", "sample", "fraction")
mac_prop_by_sample$genotype[str_detect(mac_prop_by_sample$sample, "OE")] <- "OE"
mac_prop_by_sample$genotype[str_detect(mac_prop_by_sample$sample, "WT")] <- "WT"
write.csv(mac_prop_by_sample,
          file = "results/mac-differential-abundance/mac_prop_by_sample.csv",
          row.names = F)

# miloR----
# Convert to SingleCellExperiment
mac_sce <- as.SingleCellExperiment(mac, assay = "RNA")

# Convert to milo object
mac_milo <- Milo(mac_sce)
head(colData(mac_milo))

# Build kNN graph, calculate neighborhood counts
mac_milo <- buildGraph(mac_milo, k = 40, d = 25, reduced.dim = "PCA")
mac_milo <- makeNhoods(mac_milo, prop = 0.1, k = 20, d = 25, refined = TRUE)

# Plot neighborhood size, peak should be between 50-100
p <- plotNhoodSizeHist(mac_milo)
pdf(file = "results/mac-differential-abundance/mac_milo_Nhood_size.pdf",
    useDingbats = F)
print(p)
dev.off()

# Count cells in neighborhoods
mac_milo <- countCells(mac_milo,
                       meta.data = data.frame(colData(mac_milo)),
                       sample = "sample")

# Compute neighborhood connectivity
mac_milo <- calcNhoodDistance(mac_milo, d = 25)

# Define experimental design
exp_design <- data.frame(colData(mac_milo))[, c("sample", "genotype")]
exp_design <- distinct(exp_design)
rownames(exp_design) <- exp_design$sample
exp_design

# Testing
mac_milo_res <- testNhoods(mac_milo,
  design = ~ genotype, design.df = exp_design)
mac_milo_res %>%  arrange(SpatialFDR) %>%  head()

# Inspecting and plotting results
p1 <- ggplot(mac_milo_res, aes(PValue)) + geom_histogram(bins = 50)
p2 <- ggplot(mac_milo_res, aes(logFC, -log10(SpatialFDR))) +
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)
pdf(file = "results/mac-differential-abundance/mac_milo_pvalue_spatialfdr.pdf",
    useDingbats = F)
print(p1 + p2)
dev.off()

# Build Nhood graph
mac_milo <- buildNhoodGraph(mac_milo)

# Plot DA results next to UMAP
p1 <- plotReducedDim(mac_milo,
                     dimred = "UMAP",
                     colour_by = "mac_annotation",
                     text_by = "mac_annotation",
                     text_size = 3) +
  scale_color_manual(values = colors) +
  NoLegend()
p2 <- plotNhoodGraphDA(mac_milo, mac_milo_res, alpha = 0.05)
                    # alpha here is statistical sig threshold, not transparency
pdf(file = "results/mac-differential-abundance/mac_milo_UMAP_NhoodGraph.pdf",
    height = 6.5,
    width = 12,
    useDingbats = F
    )
print(p1 + p2 + plot_layout(guides = "collect"))
dev.off()

# Annotate the Nhoods based on basic_annotation
mac_milo_res <- annotateNhoods(mac_milo,
                              mac_milo_res,
                              coldata_col = "mac_annotation")
unique(mac_milo_res$mac_annotation)
ggplot(mac_milo_res, aes(mac_annotation_fraction)) + geom_histogram(bins = 50)
mac_milo_res$mac_annotation <- ifelse(mac_milo_res$mac_annotation_fraction < 0.6,
                                    "Mixed",
                                    mac_milo_res$mac_annotation)
mac_milo_res$mac_annotation <- factor(mac_milo_res$mac_annotation,
                                    levels = maclevels)

# plot DAbeeswarm...ignored because none are significant
plotDAbeeswarm(mac_milo_res,
               group.by = "mac_annotation",
               alpha = 0.05) +
  ggtitle("Cluster assignment of DA neighborhoods") +
  theme(axis.title.y = element_blank(),
        plot.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

# Save DA results and objects
write.csv(mac_milo_res,
          file = "results/mac-differential-abundance/mac_milo_res.csv",
          row.names = F)

# Save milo object
save(mac_milo, file = "results/objects/mac_milo.Rdata")


# OCSA-DA Analysis with edgeR----
#http://bioconductor.org/books/3.14/OSCA.multisample/differential-abundance.html

# Abundances
abundances <- table(mac_sce$mac_annotation, mac_sce$sample)
abundances <- unclass(abundances)
head(abundances)

# Attaching some column metadata and making DGEList object
extra.info <- colData(mac_sce)[match(colnames(abundances), mac_sce$sample), ]
y.ab <- DGEList(abundances, samples = extra.info)
y.ab

# Filter low abundance
keep <- filterByExpr(y.ab, group = y.ab$samples$genotype)
y.ab <- y.ab[keep, ]
summary(keep)

# Design matrix
design <- model.matrix(~genotype, y.ab$samples)

# Estimate dispersion
y.ab <- estimateDisp(y.ab, design, trend = "none")
summary(y.ab$common.dispersion)
plotBCV(y.ab, cex = 1)

# QL fit
fit.ab <- glmQLFit(y.ab, design, robust = TRUE, abundance.trend = FALSE)
summary(fit.ab$var.prior)
plotQLDisp(fit.ab, cex = 1)

# Test
mac_ocsa_da_res <- glmQLFTest(fit.ab, coef = ncol(design))
summary(decideTests(mac_ocsa_da_res))
topTags(mac_ocsa_da_res, n = 34)
ocsa_da_res$table

# Workflow with normalization (assuming most labels do not change in abundance)
y.ab2 <- calcNormFactors(y.ab)
y.ab2$samples$norm.factors
y.ab2 <- estimateDisp(y.ab2, design, trend = "none")
fit.ab2 <- glmQLFit(y.ab2, design, robust = TRUE, abundance.trend = FALSE)
mac_ocsa_da_res_norm <- glmQLFTest(fit.ab2, coef = ncol(design))
topTags(mac_ocsa_da_res_norm)

# Save results
write.csv(topTags(mac_ocsa_da_res, n = 16),
          file = "results/mac-differential-abundance/mac_ocsa_da_res.csv")
write.csv(topTags(mac_ocsa_da_res_norm, n = 16),
          file = "results/mac-differential-abundance/mac_ocsa_da_res_norm.csv")

# Save sce object
save(mac_sce, file = "results/objects/mac_sce.Rdata")
