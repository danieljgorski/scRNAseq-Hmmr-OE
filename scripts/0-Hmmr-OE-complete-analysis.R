# Complete analysis of scRNAseq-Hmmr-OE data

# Make directories
output_dirs <- c(
  "results",
  "results/objects",
  "results/preprocessing",
  "results/doublet-removal",
  "results/post-clustering-qc"
)

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
source("scripts/4-Hmmr-OE-initial-clustering.R")
source("scripts/5-Hmmr-OE-post-clustering-qc.R")
#source("scripts/6-Hmmr-OE-final-clustering.R")
