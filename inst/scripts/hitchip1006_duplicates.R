# -----------------------------------------------------------------------------
# hitchip1006_duplicates.R
#
# Removes duplicated samples from `hitchip1006` (see `?miaTime::hitchip1006`).
#
# In 11 pairs of samples, the counts are identical; each pair is one subject at
# time 0 and time 0.2 with different `project` labels. The pairs are the same
# in `microbiome::atlas1006`. The time-0.2 sample of each pair is removed.
#
# Run this script from the package root.
# -----------------------------------------------------------------------------

library(SummarizedExperiment)

load("data/hitchip1006.rda")
tse <- hitchip1006

counts <- assay(tse, "counts")
key <- apply(counts, 2L, paste, collapse = ",")
pairs <- split(colnames(tse), key)
pairs <- pairs[lengths(pairs) > 1L]
stopifnot(length(pairs) == 11L, all(lengths(pairs) == 2L))
drop <- vapply(pairs, function(s) s[which.max(tse[, s]$time)], character(1))
stopifnot(all(tse[, drop]$time == 0.2))

hitchip1006 <- tse[, !colnames(tse) %in% drop]
save(hitchip1006, file = "data/hitchip1006.rda", compress = "xz")
