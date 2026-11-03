test_that("getVariogram input checks", {
    data(hitchip1006, package = "miaTime")
    tse <- hitchip1006
    expect_error(getVariogram(tse, time.col = "x", group = "subject"))
    expect_error(getVariogram(tse, time.col = "sex",
        group = "subject"))
    expect_error(getVariogram(tse, time.col = "time", group = "x"))
    expect_error(getVariogram(tse, time.col = "time",
        group = "subject", assay.type = "x"))
    expect_error(getVariogram(tse, time.col = "time",
        group = "subject", time.lag = -1))
    expect_error(getVariogram(tse, time.col = "time",
        group = "subject", time.lag = 1, window = 1))
    expect_error(getVariogram(tse, time.col = "time",
        group = "subject", time.lag = 1, time.point = "mean"))
    expect_error(getVariogram(tse, time.col = "time",
        group = "subject", time.lag = c(1, 2), time.point = c(1, 2, 3)))
    expect_error(getVariogram(tse, time.col = "time",
        group = "subject", time.lag = 1, conf.level = 1))
    # split.by must be constant within subjects
    expect_error(getVariogram(tse, time.col = "time",
        group = "subject", time.lag = 1, split.by = "time"))
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
    # At least twice as many subjects as coefficients
    pairs <- .make_quadratic_pairs(n = 5)
    res <- miaTime:::.local_variogram(pairs, lags = 10, window = 2,
        conf.level = 0.95, a0 = NA_real_)
    expect_false(res$supported)
    pairs <- .make_quadratic_pairs(n = 6)
    res <- miaTime:::.local_variogram(pairs, lags = 10, window = 2,
        conf.level = 0.95, a0 = NA_real_)
    expect_true(res$supported)
    # No pairs above L: extrapolation is not reported
    pairs <- .make_quadratic_pairs()
    pairs <- pairs[pairs$h < 10, ]
    res <- miaTime:::.local_variogram(pairs, lags = 10, window = 2,
        conf.level = 0.95, a0 = NA_real_)
    expect_false(res$supported)
})

test_that("getVariogram output", {
    data(hitchip1006, package = "miaTime")
    tse <- transformAssay(hitchip1006, method = "hellinger")
    res <- getVariogram(tse, time.col = "time", group = "subject",
        assay.type = "hellinger", method = "euclidean", squared = TRUE,
        time.lag = c(1, 3, 6))
    expect_s4_class(res, "DataFrame")
    expect_equal(res$lag, c(1, 3, 6))
    expect_true(all(res$type == "dissimilarity"))
    ok <- res$supported
    expect_true(any(ok))
    expect_true(all(res$lower[ok] < res$estimate[ok]))
    expect_true(all(res$estimate[ok] < res$upper[ok]))
    expect_true(all(is.na(res$estimate[!ok])))
    # Squaring changes the dissimilarity
    res2 <- getVariogram(tse, time.col = "time", group = "subject",
        assay.type = "hellinger", method = "euclidean", time.lag = c(1, 3, 6))
    expect_false(isTRUE(all.equal(res$estimate[ok], res2$estimate[ok])))
    # Wider intervals at a higher confidence level
    res3 <- getVariogram(tse, time.col = "time", group = "subject",
        assay.type = "hellinger", method = "euclidean", squared = TRUE,
        time.lag = c(1, 3, 6), conf.level = 0.99)
    expect_equal(res3$estimate, res$estimate)
    expect_true(all(res3$upper[ok] > res$upper[ok]))
    # Default lags
    res4 <- getVariogram(tse, time.col = "time", group = "subject",
        assay.type = "hellinger", method = "euclidean", squared = TRUE)
    expect_true(nrow(res4) > 1L)
    expect_true(all(diff(res4$lag) > 0))
})

test_that("getVariogram compares groups", {
    data(hitchip1006, package = "miaTime")
    tse <- transformAssay(hitchip1006, method = "hellinger")
    expect_message(res <- getVariogram(tse, time.col = "time",
        group = "subject", assay.type = "hellinger", method = "euclidean",
        squared = TRUE, time.lag = c(3, 6), split.by = "sex"))
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
    expect_message(res <- getVariogram(tse, time.col = "time",
        group = "subject", assay.type = "hellinger", method = "euclidean",
        squared = TRUE, time.lag = c(3, 6), split.by = "sex",
        time.point = "median"))
    expect_true("time.point" %in% colnames(res))
    expect_equal(length(unique(res$time.point[res$lag == 6])), 1L)
})

test_that("addVariogram stores results in metadata", {
    data(hitchip1006, package = "miaTime")
    tse <- transformAssay(hitchip1006, method = "hellinger")
    args <- list(time.col = "time", group = "subject",
        assay.type = "hellinger", method = "euclidean", time.lag = c(3, 6))
    res <- do.call(getVariogram, c(list(x = tse), args))
    tse <- do.call(addVariogram, c(list(x = tse), args))
    expect_equal(metadata(tse)[["variogram"]], res)
    tse <- do.call(addVariogram,
        c(list(x = tse, name = "vario"), args))
    expect_equal(metadata(tse)[["vario"]], res)
})

test_that("Permutation p-values of comparisons", {
    data(hitchip1006, package = "miaTime")
    tse <- transformAssay(hitchip1006, method = "hellinger")
    args <- list(tse, time.col = "time", group = "subject",
        assay.type = "hellinger", method = "euclidean", squared = TRUE,
        time.lag = c(1, 3), split.by = "sex")
    expect_error(do.call(getVariogram, c(args, permutations = -1)))
    expect_error(do.call(getVariogram, c(args, permutations = 1.5)))
    # No permutation column by default
    res <- do.call(getVariogram, args)
    expect_false("p.value.permutation" %in% colnames(res))
    set.seed(1)
    res <- do.call(getVariogram, c(args, permutations = 19))
    expect_true("p.value.permutation" %in% colnames(res))
    # Only comparisons get a permutation p-value
    expect_true(all(is.na(
        res$p.value.permutation[res$type == "dissimilarity"])))
    p <- res$p.value.permutation[res$type == "ratio" & res$supported]
    expect_true(length(p) > 0L)
    expect_true(all(p >= 1 / 20 & p <= 1))
    # The estimates do not change
    res0 <- do.call(getVariogram, args)
    expect_equal(res$estimate, res0$estimate)
    expect_equal(res$p.value, res0$p.value)
    # Reproducible with a seed
    set.seed(1)
    res2 <- do.call(getVariogram, c(args, permutations = 19))
    expect_equal(res$p.value.permutation, res2$p.value.permutation)
})
