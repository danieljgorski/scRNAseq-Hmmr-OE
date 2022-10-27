# Re-clustering cleaned object, using reference based integration with RPCA

# Load libraries
library(Seurat) #v4.0.1

# Load object
load("results/objects/obj_clean.Rdata")

# Split object on samples
obj_list <- SplitObject(obj, split.by = "sample")

# SCTransform objects individually with sped up glmGamPoi method
obj_list <- lapply(X = obj_list, FUN = SCTransform, method = "glmGamPoi")

# Find integration features, remove Hmmr and FIJ5 (manual over-expression)
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)
features <- features[features != "Hmmr"]
features <- features[features != "FIJ5"]

# Prep for SCT integration, run PCA on individual samples
obj_list <- PrepSCTIntegration(object.list = obj_list,
                               anchor.features = features)
obj_list <- lapply(X = obj_list, FUN = RunPCA, features = features)

# Find integration anchors, using WT samples as reference
controls <- c(2, 3, 6, 10, 11, 14) # indices of WT samples
obj_anchors <- FindIntegrationAnchors(object.list = obj_list,
                                      normalization.method = "SCT",
                                      reference = controls,
                                      anchor.features = features,
                                      dims = 1:25,
                                      reduction = "rpca",
                                      k.anchor = 5)

# Save anchors, due to the large memory requirements of integration
save(obj_anchors, file = "results/objects/obj_anchors_clean.Rdata")
remove(obj)
remove(obj_list)

# Integrate data
obj <- IntegrateData(anchorset = obj_anchors,
                     normalization.method = "SCT",
                     dims = 1:25)

# Run PCA, UMAP and cluster integrated object
obj <- RunPCA(obj, npcs = 50, verbose = T)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:25, verbose = T)
obj <- FindNeighbors(obj, dims = 1:25, verbose = T)
obj <- FindClusters(obj, resolution = 0.4, verbose = T)
DimPlot(obj, label = T)

# Normalizing and scaling RNA assay
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj,
                     normalization.method = "LogNormalize",
                     scale.factor = 10000,
                     verbose = T)
obj <- ScaleData(obj, features = rownames(obj), verbose = T)

# Factor genotype level
obj@meta.data$genotype <- factor(obj@meta.data$genotype, levels = c("WT", "OE"))

# Save object
save(obj, file = "results/objects/obj_integrated_clean.Rdata")

# Clear memory
rm(list = ls())
gc()
