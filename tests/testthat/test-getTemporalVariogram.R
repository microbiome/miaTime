test_that("getTemporalVariogram input checks", {
    data(hitchip1006, package = "miaTime")
    tse <- hitchip1006
    expect_error(getTemporalVariogram(tse, time.col = "x", group = "subject"))
    expect_error(getTemporalVariogram(tse, time.col = "sex",
        group = "subject"))
    expect_error(getTemporalVariogram(tse, time.col = "time", group = "x"))
    expect_error(getTemporalVariogram(tse, time.col = "time",
        group = "subject", assay.type = "x"))
    expect_error(getTemporalVariogram(tse, time.col = "time",
        group = "subject", lags = -1))
    expect_error(getTemporalVariogram(tse, time.col = "time",
        group = "subject", lags = 1, window = 1))
    expect_error(getTemporalVariogram(tse, time.col = "time",
        group = "subject", lags = 1, at.time = "mean"))
    expect_error(getTemporalVariogram(tse, time.col = "time",
        group = "subject", lags = c(1, 2), at.time = c(1, 2, 3)))
    expect_error(getTemporalVariogram(tse, time.col = "time",
        group = "subject", lags = 1, conf.level = 1))
    # split.by must be constant within subjects
    expect_error(getTemporalVariogram(tse, time.col = "time",
        group = "subject", lags = 1, split.by = "time"))
})

# Pairs of n subjects with lags spread around 10 and d an exact quadratic in
# log(h / 10)
.make_quadratic_pairs <- function(n = 20, at = 0.3){
    h <- rep(c(6, 8, 10, 12, 15), n)
    u <- log(h / 10)
    data.frame(
        subject = rep(paste0("s", seq_len(n)), each = 5),
        h = h, mid = 10, d = at + 0.1 * u - 0.05 * u^2)
}

test_that("Local estimate is exact when the variogram is quadratic", {
    pairs <- .make_quadratic_pairs()
    res <- miaTime:::.local_variogram(pairs, lags = 10, window = 2,
        conf.level = 0.95, a0 = NA_real_)
    expect_equal(res$estimate, 0.3)
    expect_true(res$supported)
    expect_equal(res$n.subjects, 20L)
    expect_equal(res$df, 19)
    expect_equal(res$se, 0, tolerance = 1e-8)
})

test_that("D(L) is reported only at supported lags", {
    # Fewer than five subjects
    pairs <- .make_quadratic_pairs(n = 4)
    res <- miaTime:::.local_variogram(pairs, lags = 10, window = 2,
        conf.level = 0.95, a0 = NA_real_)
    expect_false(res$supported)
    expect_true(is.na(res$estimate))
    # No pairs above L: extrapolation is not reported
    pairs <- .make_quadratic_pairs()
    pairs <- pairs[pairs$h < 10, ]
    res <- miaTime:::.local_variogram(pairs, lags = 10, window = 2,
        conf.level = 0.95, a0 = NA_real_)
    expect_false(res$supported)
})

test_that("getTemporalVariogram output", {
    data(hitchip1006, package = "miaTime")
    tse <- transformAssay(hitchip1006, method = "hellinger")
    res <- getTemporalVariogram(tse, time.col = "time", group = "subject",
        assay.type = "hellinger", method = "euclidean", squared = TRUE,
        lags = c(1, 3, 6))
    expect_s4_class(res, "DataFrame")
    expect_equal(res$lag, c(1, 3, 6))
    expect_true(all(res$type == "dissimilarity"))
    ok <- res$supported
    expect_true(any(ok))
    expect_true(all(res$lower[ok] < res$estimate[ok]))
    expect_true(all(res$estimate[ok] < res$upper[ok]))
    expect_true(all(is.na(res$estimate[!ok])))
    # Squaring changes the dissimilarity
    res2 <- getTemporalVariogram(tse, time.col = "time", group = "subject",
        assay.type = "hellinger", method = "euclidean", lags = c(1, 3, 6))
    expect_false(isTRUE(all.equal(res$estimate[ok], res2$estimate[ok])))
    # Wider intervals at a higher confidence level
    res3 <- getTemporalVariogram(tse, time.col = "time", group = "subject",
        assay.type = "hellinger", method = "euclidean", squared = TRUE,
        lags = c(1, 3, 6), conf.level = 0.99)
    expect_equal(res3$estimate, res$estimate)
    expect_true(all(res3$upper[ok] > res$upper[ok]))
    # Default lags
    res4 <- getTemporalVariogram(tse, time.col = "time", group = "subject",
        assay.type = "hellinger", method = "euclidean", squared = TRUE)
    expect_true(nrow(res4) > 1L)
    expect_true(all(diff(res4$lag) > 0))
})

test_that("getTemporalVariogram compares groups", {
    data(hitchip1006, package = "miaTime")
    tse <- transformAssay(hitchip1006, method = "hellinger")
    expect_message(res <- getTemporalVariogram(tse, time.col = "time",
        group = "subject", assay.type = "hellinger", method = "euclidean",
        squared = TRUE, lags = c(3, 6), split.by = "sex"))
    expect_true("sex" %in% colnames(res))
    expect_equal(nrow(res), 6L)
    est <- res[res$type == "dissimilarity", ]
    rat <- res[res$type == "ratio", ]
    expect_equal(unique(rat$sex), "male / female")
    for( l in c(3, 6) ){
        e <- est[est$lag == l, ]
        r <- rat[rat$lag == l, ]
        expect_equal(r$estimate,
            e$estimate[e$sex == "male"] / e$estimate[e$sex == "female"])
        expect_equal(r$n.subjects, sum(e$n.subjects))
        expect_true(r$p.value >= 0 && r$p.value <= 1)
    }
    expect_true(all(is.na(est$p.value)))
    # Comparison at a common time
    expect_message(res <- getTemporalVariogram(tse, time.col = "time",
        group = "subject", assay.type = "hellinger", method = "euclidean",
        squared = TRUE, lags = c(3, 6), split.by = "sex",
        at.time = "median"))
    expect_true("at.time" %in% colnames(res))
    expect_equal(length(unique(res$at.time[res$lag == 6])), 1L)
})

test_that("addTemporalVariogram stores results in metadata", {
    data(hitchip1006, package = "miaTime")
    tse <- transformAssay(hitchip1006, method = "hellinger")
    args <- list(time.col = "time", group = "subject",
        assay.type = "hellinger", method = "euclidean", lags = c(3, 6))
    res <- do.call(getTemporalVariogram, c(list(x = tse), args))
    tse <- do.call(addTemporalVariogram, c(list(x = tse), args))
    expect_equal(metadata(tse)[["temporal_variogram"]], res)
    tse <- do.call(addTemporalVariogram,
        c(list(x = tse, name = "vario"), args))
    expect_equal(metadata(tse)[["vario"]], res)
})
