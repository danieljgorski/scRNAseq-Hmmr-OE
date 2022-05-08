# Preprocessing count matrices

# Load libraries
library(Seurat) #v4.0.1
library(dplyr)

# Read in 10x outputs and create Seurat objects
day_3 <- Read10X("data/count_matrices/day_3_filtered_feature_bc_matrix")
day_3 <- CreateSeuratObject(counts = day_3,
                            project = "474",
                            min.cells = 3,
                            min.features = 200,
                            names.field = 2,
                            names.delim = "-")
day_7 <- Read10X("data/count_matrices/day_7_filtered_feature_bc_matrix")
day_7 <- CreateSeuratObject(counts = day_7,
                            project = "153",
                            min.cells = 3,
                            min.features = 200,
                            names.field = 2,
                            names.delim = "-")

# Add metadata
day_3 <- AddMetaData(day_3, metadata = "day_3", col.name = "timepoint")
day_3@meta.data <- day_3@meta.data %>%
  mutate(sample =
           case_when(orig.ident == "1" ~ "S47_B1_OE",
                     orig.ident == "2" ~ "S47_B3_OE",
                     orig.ident == "3" ~ "S47_B5_WT",
                     orig.ident == "4" ~ "S47_B6_WT",
                     orig.ident == "5" ~ "S48_R1_OE",
                     orig.ident == "6" ~ "S48_R3_OE",
                     orig.ident == "7" ~ "S48_R5_WT"))
day_3@meta.data <- day_3@meta.data %>%
  mutate(date =
           case_when(orig.ident == "1" ~ "16_Aug_21_day_3",
                     orig.ident == "2" ~ "16_Aug_21_day_3",
                     orig.ident == "3" ~ "16_Aug_21_day_3",
                     orig.ident == "4" ~ "16_Aug_21_day_3",
                     orig.ident == "5" ~ "19_Aug_21_day_3",
                     orig.ident == "6" ~ "19_Aug_21_day_3",
                     orig.ident == "7" ~ "19_Aug_21_day_3"))
day_3@meta.data <- day_3@meta.data %>%
  mutate(genotype =
           case_when(orig.ident == "1" ~ "OE",
                     orig.ident == "2" ~ "OE",
                     orig.ident == "3" ~ "WT",
                     orig.ident == "4" ~ "WT",
                     orig.ident == "5" ~ "OE",
                     orig.ident == "6" ~ "OE",
                     orig.ident == "7" ~ "WT"))
day_7 <- AddMetaData(day_7, metadata = "day_7", col.name = "timepoint")
day_7@meta.data <- day_7@meta.data %>%
  mutate(sample =
           case_when(orig.ident == "1" ~ "S18_G2_WT",
                     orig.ident == "2" ~ "S18_G1_OE",
                     orig.ident == "3" ~ "S18_G3_WT",
                     orig.ident == "4" ~ "S18_R1_OE",
                     orig.ident == "5" ~ "S18_G8_WT",
                     orig.ident == "6" ~ "S18_G6_OE",
                     orig.ident == "7" ~ "S18_G7_OE"))
day_7 <- AddMetaData(day_7, metadata = "22_Jun_20_day_7", col.name = "date")
day_7@meta.data <- day_7@meta.data %>%
  mutate(genotype =
           case_when(orig.ident == "1" ~ "WT",
                     orig.ident == "2" ~ "OE",
                     orig.ident == "3" ~ "WT",
                     orig.ident == "4" ~ "OE",
                     orig.ident == "5" ~ "WT",
                     orig.ident == "6" ~ "OE",
                     orig.ident == "7" ~ "OE"))

# Merge objects and remove singles
obj <- merge(x = day_3, y = day_7)
remove(day_3)
remove(day_7)

# Quality control filtering
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
cells_pre_qc <- length(colnames(obj))

# Pre-filter QC metrics by sample
p <- VlnPlot(obj,
             features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
             group.by = "sample",
             pt.size = 0)
pdf(file = "results/preprocessing/VlnPlot_QC_metrics_pre-filter.pdf",
    width = 12,
    height = 6,
    pointsize = 12,
    useDingbats = F)
print(p)
dev.off()

# Filter
obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 10)
cells_post_qc <- length(colnames(obj))
percent_passed <- (cells_post_qc / cells_pre_qc) * 100
qc_filter <- as.data.frame(cells_pre_qc)
qc_filter$cells_post_qc <- cells_post_qc
qc_filter$percent_passed <- percent_passed
write.csv(qc_filter,
          file = "results/preprocessing/qc_filter.csv",
          row.names = F)

# Post-filter QC metrics by sample
p <- VlnPlot(obj,
             features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
             group.by = "sample",
             pt.size = 0)
pdf(file = "results/preprocessing/VlnPlot_QC_metrics_post-filter.pdf",
    width = 12,
    height = 6,
    pointsize = 12,
    useDingbats = F)
print(p)
dev.off()

# Save object
save(obj, file = "results/objects/obj.Rdata")
