# Re-clustering macrophages using using reference based integration with RPCA

# Load libraries
library(Seurat) #v4.0.1

# Load full object
load("results/objects/mac.Rdata")

# Split object on samples
obj_list <- SplitObject(mac, split.by = "sample")

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
mac_anchors <- FindIntegrationAnchors(object.list = obj_list,
                                      normalization.method = "SCT",
                                      reference = controls,
                                      anchor.features = features,
                                      dims = 1:25,
                                      reduction = "rpca",
                                      k.anchor = 5)

# Save anchors, due to the large memory requirements of integration
save(mac_anchors, file = "results/objects/mac_anchors_clean.Rdata")
remove(mac)
remove(obj_list)

# Integrate data
mac <- IntegrateData(anchorset = mac_anchors,
                     normalization.method = "SCT",
                     dims = 1:25)

# Run PCA, UMAP and cluster integrated object
mac <- RunPCA(mac, npcs = 50, verbose = T)
ElbowPlot(mac, ndims = 50)
mac <- RunUMAP(mac, reduction = "pca", dims = 1:25, verbose = T)
mac <- FindNeighbors(mac, dims = 1:25, verbose = T)
mac <- FindClusters(mac, resolution = 0.6, verbose = T) # reduced resolution
DimPlot(mac, label = T)

# Normalizing and scaling RNA assay
DefaultAssay(mac) <- "RNA"
mac <- NormalizeData(mac,
                     normalization.method = "LogNormalize",
                     scale.factor = 10000,
                     verbose = T)
mac <- ScaleData(mac, features = rownames(mac), verbose = T)

# Factor genotype level
mac@meta.data$genotype <- factor(mac@meta.data$genotype, levels = c("WT", "OE"))

# Save object
save(mac, file = "results/objects/mac_integrated_clean.Rdata")
