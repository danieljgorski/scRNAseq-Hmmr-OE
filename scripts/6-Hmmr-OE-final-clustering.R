# Re-clustering cleaned object

# Load libraries
library(Seurat) #v4.0.1

# Load object
load("results/objects/obj.Rdata")

# Exclude low-quality clusters determined during post-clustering-qc
obj <- subset(x = obj, idents = "33", invert = TRUE)

# Remove extra meta data columns
obj@meta.data <- select(obj@meta.data,
                        -starts_with("pANN"),
                        -c(nCount_SCT,
                           nFeature_SCT,
                           SCT_snn_res.0.8,
                           seurat_clusters,
                           doublet_classification))

# Keep only necessary counts and data slots
obj <- DietSeurat(obj,
                  counts = T,
                  data = T,
                  scale.data = F,
                  features = NULL,
                  assays = "RNA",
                  dimreducs = NULL,
                  graphs = NULL)

# Save pruned object
save(obj, file = "results/objects/obj.Rdata")

# Clear memory
rm(list = ls())
gc()

# Load cleaned object
load("results/objects/obj.Rdata")

# SCTransform
obj <- SCTransform(obj, verbose = T, variable.features.n = 3002)

# Take Hmmr and FIJ5 out of variable features
var_features <- VariableFeatures(obj, assay = "SCT")
var_features <- var_features[var_features != "Hmmr"]
var_features <- var_features[var_features != "FIJ5"]
VariableFeatures(obj, assay = "SCT") <- var_features

# PCA
obj <- RunPCA(obj, npcs = 50, verbose = T)
ElbowPlot(obj, ndims = 50)

# UMAP and clustering
DefaultAssay(obj) <- "SCT"
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, verbose = T)
obj <- FindNeighbors(obj, dims = 1:30, verbose = T)
obj <- FindClusters(obj, resolution = 0.8, verbose = T)
DimPlot(obj, label = T, raster = F)

# Normalizing and scaling RNA assay
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj,
                     normalization.method = "LogNormalize",
                     scale.factor = 10000,
                     verbose = T)
obj <- ScaleData(obj, features = rownames(obj), verbose = T)

# Save object
save(obj, file = "results/objects/obj.Rdata")
