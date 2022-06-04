# Complete analysis of scRNAseq-Hmmr-OE data

# Make directories
output_dirs <- c("results",
                 "results/basic-annotation",
                 "results/cluster-markers",
                 "results/differential-abundance",
                 "results/differential-gene-expression",
                 "results/dimplots",
                 "results/doublet-removal",
                 "results/mac-annotation",
                 "results/mac-cluster-markers",
                 "results/mac-common-dge",
                 "results/mac-differential-abundance",
                 "results/mac-differential-gene-expression",
                 "results/mac-dimplots",
                 "results/mac-ora",
                 "results/mac-violin-plots",
                 "results/mac-volcano-plots",
                 "results/objects",
                 "results/ora",
                 "results/post-clustering-qc",
                 "results/preprocessing",
                 "results/volcano-plots")

for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}

# Run analysis in order
source("scripts/1-Hmmr-OE-preprocessing.R")
source("scripts/2-Hmmr-OE-doublet-classification.R")
source("scripts/3-Hmmr-OE-doublet-removal.R")
source("scripts/4-Hmmr-OE-initial-integrated-clustering.R")
source("scripts/5-Hmmr-OE-post-clustering-qc.R")
source("scripts/6-Hmmr-OE-final-clustering.R")
source("scripts/7-Hmmr-OE-basic-annotation.R")
source("scripts/8-Hmmr-OE-dimplots.R")
source("scripts/9-Hmmr-OE-cluster-markers.R")
source("scripts/10-Hmmr-OE-differential-abundance.R")
source("scripts/11-Hmmr-OE-differential-gene-expression.R")
source("scripts/12-Hmmr-OE-volcano-plots.R")
source("scripts/13-Hmmr-OE-ora.R")
source("scripts/14-Hmmr-OE-mac-subset.R")
source("scripts/15-Hmmr-OE-mac-integrated-clustering.R")
source("scripts/16-Hmmr-OE-mac-annotation.R")
source("scripts/17-Hmmr-OE-mac-dimplots.R")
source("scripts/18-Hmmr-OE-mac-cluster-markers.R")
source("scripts/19-Hmmr-OE-mac-differential-abundance.R")
source("scripts/20-Hmmr-OE-mac-differential-gene-expression.R")
source("scripts/21-Hmmr-OE-mac-volcano-plots.R")
source("scripts/22-Hmmr-OE-mac-ora.R")
source("scripts/23-Hmmr-OE-mac-common-dge.R")
source("scripts/24-Hmmr-OE-mac-violin-plots.R")
