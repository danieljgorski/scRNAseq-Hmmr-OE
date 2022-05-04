# make directories
output_dirs <- c(
  "results",
  "results/objects"
)

for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}

# run analysis in order
source("scripts/1-Hmmr-OE-pre-processing.R")
