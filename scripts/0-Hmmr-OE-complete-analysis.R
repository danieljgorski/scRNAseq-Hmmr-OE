# Complete analysis of scRNAseq-Hmmr-OE data

# Run analysis in order
source("scripts/1-Hmmr-OE-preprocessing.R")
source("scripts/2-Hmmr-OE-doublet-classification.R")
source("scripts/3-Hmmr-OE-doublet-removal.R")
source("scripts/4-Hmmr-OE-initial-integrated-clustering.R")
#source("scripts/5-Hmmr-OE-post-clustering-qc.R")
#source("scripts/6-Hmmr-OE-final-clustering.R")
#source("scripts/6-Hmmr-OE-basic-annotation.R")
