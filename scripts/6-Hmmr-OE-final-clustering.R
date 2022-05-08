# Load libraries
library(Seurat) #v4.0.1

# Load object
load("results/objects/obj.Rdata")

# Re-clustering cleaned object

# SCTransform
obj <- SCTransform(obj, verbose = T, variable.features.n = 3002)

#Taking Hmmr and FIJ5 out of variable features
var.features <- VariableFeatures(obj, assay = "SCT")
which(var.features=="Hmmr")
which(var.features=="FIJ5")
var.features <- var.features[var.features!= "Hmmr"]
var.features <- var.features[var.features!= "FIJ5"]
which(var.features=="Hmmr")
which(var.features=="FIJ5")
VariableFeatures(obj, assay = "SCT") <- var.features

# PCA
obj <- RunPCA(obj, npcs = 50, verbose = T)
ElbowPlot(obj, ndims = 50)

# UMAP and clustering
DefaultAssay(obj) <- "SCT"
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, verbose = T)
obj <- FindNeighbors(obj, dims = 1:30,verbose = T)
obj <- FindClusters(obj, resolution = 0.8, verbose = T)
DimPlot(obj, label = T, raster = F)

# Normalizing and scaling RNA assay
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj,
                     normalization.method = "LogNormalize",
                     scale.factor = 10000,
                     verbose = T)
obj <- ScaleData(obj, features = rownames(obj), verbose = T)

## Save object
save(obj, file = "results/objects/obj.Rdata")
