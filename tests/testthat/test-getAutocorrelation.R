# A series with known autocorrelation, for checks real data cannot give.
set.seed(1)
n <- 300
ar <- as.numeric(arima.sim(list(ar = 0.7), n))
se <- SummarizedExperiment(
    assays = list(counts = rbind(ar = ar)),
    colData = DataFrame(t = seq_len(n) * 2))

# Artificial gut vessels sampled every hour; hours from the sample names.
data(SilvermanAGutData)
tse <- SilvermanAGutData
tse <- tse[, tse$SampleType == "Hourly" &
    tse$Normal_Noise_Sample == "Normal"]
tse <- tse[, grepl("^T\\.+[0-9]+d\\.[0-9]+h", colnames(tse)) &
    colSums(assay(tse, "counts")) > 0]
days <- as.numeric(sub("^T\\.+([0-9]+)d.*", "\\1", colnames(tse)))
hours <- as.numeric(sub("^T\\.+[0-9]+d\\.([0-9]+)h.*", "\\1", colnames(tse)))
tse$hour <- days * 24 + hours
tse <- tse[1:20, ]

test_that("getAutocorrelation checks its input", {
    expect_error(getAutocorrelation(se, time.col = "missing"))
    expect_error(getAutocorrelation(tse, time.col = "SampleType"), "numeric")
    expect_error(getAutocorrelation(se, time.col = "t", max.lag = "a"))
    irregular <- se[, 1:3]
    irregular$t <- c(0, 1, 2.5)
    expect_error(getAutocorrelation(irregular, time.col = "t"),
        "regular grid")
})

test_that("getAutocorrelation equals stats::acf without gaps", {
    res <- getAutocorrelation(se, time.col = "t", max.lag = 3)
    expect_equal(res[["acf"]],
        as.vector(acf(ar, lag.max = 3, plot = FALSE)[["acf"]][2:4]))
    expect_equal(res[["lag"]], 1:3)
    expect_equal(res[["time.lag"]], c(2, 4, 6))
    expect_equal(res[["n.pairs"]], n - 1:3)
    expect_false("group" %in% names(res))
})

test_that("getAutocorrelation treats missing time points as missing", {
    keep <- seq_len(n) %% 5 != 0
    res <- getAutocorrelation(se[, keep], time.col = "t", max.lag = 2)
    # Pairs exist only where both samples are present.
    expect_equal(res[["n.pairs"]], c(
        sum(keep[-n] & keep[-1]), sum(keep[seq_len(n - 2)] & keep[-(1:2)])))
    expect_equal(res[["acf"]][1L], 0.7, tolerance = 0.1)
})

test_that("getAutocorrelation works per subject and averages replicates", {
    expect_message(
        res <- getAutocorrelation(tse, time.col = "hour", group = "Vessel",
            max.lag = 5),
        "replicated samples")
    expect_true("Vessel" %in% names(res))
    expect_setequal(unique(res[["Vessel"]]), as.character(unique(tse$Vessel)))
    expect_identical(nrow(res), nrow(tse) * 4L * 5L)
    expect_true(all(res[["acf"]] >= -1 & res[["acf"]] <= 1, na.rm = TRUE))
})

test_that("addAutocorrelation stores the table in metadata", {
    x <- suppressMessages(addAutocorrelation(tse, time.col = "hour",
        group = "Vessel", max.lag = 3, name = "acf_vessel"))
    expect_equal(metadata(x)[["acf_vessel"]],
        suppressMessages(getAutocorrelation(tse, time.col = "hour",
            group = "Vessel", max.lag = 3)))
})
