# Load libraries
library(Seurat) #v4.0.1

# Load object
load("results/objects/obj.Rdata")

# Doublet removal

# Cluster merged object to visualize doublets
obj <- SCTransform(obj, verbose = T)
obj <- RunPCA(obj, verbose = T)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, verbose = T)
obj <- FindNeighbors(obj, dims = 1:30, verbose = T)
obj <- FindClusters(obj, resolution = 0.8, verbose = T)
obj@meta.data$doublet_classification <- factor(obj@meta.data$doublet_classification,
                                               levels = c("Singlet", "Doublet"))
p <- DimPlot(obj, group.by = "doublet_classification", raster = F)
pdf(file = "results/doublet-removal/DimPlot_doublet_classification.pdf",
    useDingbats = F)
print(p)
dev.off()

# Subset for singlets and save object
obj <- subset(x = obj, subset = doublet_classification == "Singlet")
save(obj, file = "results/objects/obj.Rdata")
