# Subsetting fibroblast clusters for re-clustering and targeted analysis

# Load libraries
library(Seurat) #v4.0.1
library(ggplot2)
library(patchwork)

# Load full object
load("results/objects/obj_annotated.Rdata")

# Subset fibroblast clusters
fibro <- subset(x = obj, idents = c("Fibro-1",
                                    "Fibro-2",
                                    "Fibro-3"))

# Save fibroblast object
save(fibro, file = "results/objects/fibro.Rdata")
