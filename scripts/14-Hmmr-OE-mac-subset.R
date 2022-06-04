# Subsetting macrophage clusters for re-clustering and targeted analysis

# Load libraries
library(Seurat) #v4.0.1

# Load full object
load("results/objects/obj_annotated.Rdata")

# Subset macrophage clusters
mac <- subset(x = obj, idents = c("Mac-1",
                                   "Mac-2",
                                   "Mac-3",
                                   "Mac-4",
                                   "Mac-5",
                                   "Mac-6",
                                   "Mac-7",
                                   "Mac-8",
                                   "Mac-9"))
mac <- DietSeurat(mac, assays = "RNA")

# Subset macrophages inside Cycling-1
cyc <- subset(x = obj, idents = "Cycling-1")
cyc_mac <- subset(x = cyc, subset = Fcgr1 > 0.1 & Adgre1 > 0.1 & Cd68 > 0.5)
cyc_mac <- DietSeurat(cyc_mac, assays = "RNA")

# Merge all macrophages
mac <- merge(x = mac, y = cyc_mac)

# Save macrophage object
save(mac, file = "results/objects/mac.Rdata")
