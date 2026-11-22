# -----------------------------------------------------------------------------
# SilvermanAGutData_time.R
#
# Adds the sampling time `Time.hr` (hours from midnight of day 0) to
# `SilvermanAGutData` (see `?miaTime::SilvermanAGutData`).
#
# The time is only given in the sample names:
#   - hourly samples, for example "T..19d.20h.00m...": day 19 at 20:00;
#   - daily samples, for example "T..25.Days..." or "T..0.V.2...": day 25 or
#     day 0. Daily samples were collected at 15:00 (Silverman et al. 2018,
#     Additional file 4).
# Seven samples whose description notes a wrong tube get NA, because their
# day is uncertain.
#
# Run this script from the package root.
# -----------------------------------------------------------------------------

library(SummarizedExperiment)

load("data/SilvermanAGutData.rda")
tse <- SilvermanAGutData
sample <- colnames(tse)

hourly <- grepl("^T\\.+[0-9]+d\\.[0-9]+h", sample)
daily <- grepl("^T\\.+[0-9]+\\.(Days?|V)", sample)
stopifnot(!any(hourly & daily), all(hourly | daily))

time <- rep(NA_real_, ncol(tse))
time[hourly] <-
    24 * as.numeric(sub("^T\\.+([0-9]+)d.*", "\\1", sample[hourly])) +
    as.numeric(sub("^T\\.+[0-9]+d\\.([0-9]+)h.*", "\\1", sample[hourly]))
time[daily] <-
    24 * as.numeric(sub("^T\\.+([0-9]+)\\..*", "\\1", sample[daily])) + 15
time[grepl("wrong ?tube", tse$Description, ignore.case = TRUE)] <- NA

tse$Time.hr <- time
SilvermanAGutData <- tse
save(SilvermanAGutData, file = "data/SilvermanAGutData.rda", compress = "xz")
