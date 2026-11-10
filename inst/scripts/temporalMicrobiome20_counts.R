# -----------------------------------------------------------------------------
# temporalMicrobiome20_counts.R
#
# Rebuilds the `counts` assay of `temporalMicrobiome20` (see
# `?miaTime::temporalMicrobiome20`) from the quantitative microbiome profiles
# in Supplementary Data 1, sheet S1-4, of Vandeputte et al. (2021):
# https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-021-27098-7/MediaObjects/41467_2021_27098_MOESM3_ESM.xlsx
#
# Stored as integers, values above .Machine$integer.max (about 2.1e9 cells
# per gram) became NA. The values are rounded and stored as doubles.
#
# Run this script from the package root with the downloaded file in it.
# -----------------------------------------------------------------------------

library(readxl)
library(SummarizedExperiment)

load("data/temporalMicrobiome20.rda")
tse <- temporalMicrobiome20

qmp <- as.data.frame(read_excel(
    "41467_2021_27098_MOESM3_ESM.xlsx", sheet = "S1-4", na = "NA"))
counts <- t(vapply(rownames(tse), function(genus) as.numeric(qmp[[genus]]),
    numeric(nrow(qmp))))
colnames(counts) <- paste0("Sample_", qmp[[1]])
counts <- round(counts[, colnames(tse)])

assay(tse, "counts") <- counts
temporalMicrobiome20 <- tse
save(temporalMicrobiome20, file = "data/temporalMicrobiome20.rda",
    compress = "xz")
