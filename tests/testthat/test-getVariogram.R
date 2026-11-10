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

test_that("Bounded dissimilarities are analysed on the -log(1 - d) scale", {
    data(hitchip1006, package = "miaTime")
    tse <- transformAssay(hitchip1006, method = "relabundance")
    args <- list(tse, time.col = "time", group = "subject",
        assay.type = "relabundance", method = "bray", time.lag = c(1, 3))
    pairs <- miaTime:::.get_within_subject_pairs(tse, "time", "subject",
        "relabundance", "bray", NULL, FALSE, FALSE, list())
    pairs$d <- -log(1 - pairs$d)
    ref <- miaTime:::.local_variogram(pairs, c(1, 3), 2, 0.95, c(NA, NA))
    res <- do.call(getVariogram, c(args, bounded = TRUE))
    expect_equal(res$estimate, ref$estimate)
    expect_true(any(res$supported))
    expect_true(all(res$estimate > do.call(getVariogram, args)$estimate,
        na.rm = TRUE))
    # Dissimilarities above 1 are not bounded
    expect_error(getVariogram(tse, time.col = "time", group = "subject",
        assay.type = "counts", method = "euclidean", time.lag = 1,
        bounded = TRUE))
    expect_error(do.call(getVariogram, c(args, bounded = NA)))
})

test_that("Relative change divides D(L) by the dissimilarity between subjects", {
    data(hitchip1006, package = "miaTime")
    tse <- transformAssay(hitchip1006, method = "relabundance")
    # Subjects with two or more samples, and ten with one
    n <- table(tse$subject)
    tse <- tse[, tse$subject %in% c(names(n)[n > 1], names(n)[n == 1][1:10])]
    args <- list(tse, time.col = "time", group = "subject",
        assay.type = "relabundance", method = "bray", time.lag = c(1, 3),
        split.by = "sex")
    res <- do.call(getVariogram, c(args, relative = TRUE))
    # Every pair of subjects counts equally
    d <- as.matrix(getDissimilarity(t(assay(tse, "relabundance")),
        method = "bray"))
    s <- tse$subject
    sp <- combn(unique(s), 2)
    B <- mean(apply(sp, 2, function(p) mean(d[s == p[1], s == p[2]])))
    est <- res$type == "dissimilarity"
    expect_equal(res$relative[est], res$estimate[est] / B)
    expect_true(all(is.na(res$relative[!est])))
    expect_false("relative" %in% names(do.call(getVariogram, args)))
    # On the same scale as D(L)
    res <- do.call(getVariogram, c(args, relative = TRUE, bounded = TRUE))
    B <- mean(apply(sp, 2, function(p)
        mean(-log(1 - d[s == p[1], s == p[2]]))))
    expect_equal(res$relative[est], res$estimate[est] / B)
})

test_that("Pairs are classified by a sample-level state", {
    # Every subject has pairs of both kinds; D(L) is exactly linear in
    # log(h / 10), twice as high for pairs starting in state b
    h <- rep(c(6, 8, 10, 12, 15), 2)
    u <- log(h / 10)
    pairs <- data.frame(
        subject = rep(paste0("s", 1:12), each = 10),
        h = rep(h, 12), mid = 10,
        level = factor(rep(rep(c("a", "b"), each = 5), 12)))
    pairs$d <- ifelse(pairs$level == "a", 0.3 + 0.1 * u, 0.6 + 0.05 * u)
    e <- lapply(c("a", "b"), function(g) miaTime:::.local_variogram(
        pairs[pairs$level == g, ], 10, 2, 0.95, NA_real_))
    res <- miaTime:::.compare_pair_levels(pairs, "a", "b", 10, 2, NA_real_,
        e[[1]], e[[2]])
    expect_equal(res$estimate, 2)
    expect_equal(res$se, 0, tolerance = 1e-8)
    expect_equal(res$n.subjects, 12L)
    expect_equal(res$df, 11)

    data(hitchip1006, package = "miaTime")
    tse <- transformAssay(hitchip1006, method = "hellinger")
    tse$state <- ifelse(tse$time < 2, "early", "late")
    args <- list(tse, time.col = "time", group = "subject",
        assay.type = "hellinger", method = "euclidean", squared = TRUE,
        time.lag = c(1, 3), split.by = "state")
    expect_message(res <- do.call(getVariogram, args), "earlier sample")
    expect_identical(unique(res$state),
        c("early", "late", "late / early"))
    # Each pair takes the state of its earlier sample
    pairs <- miaTime:::.get_within_subject_pairs(tse, "time", "subject",
        "hellinger", "euclidean", "state", TRUE, FALSE, list())
    expect_true(attr(pairs, "pair.level"))
    expect_true(all(pairs$level[pairs$mid - pairs$h / 2 >= 2] == "late"))
    expect_true(all(pairs$level[pairs$mid - pairs$h / 2 < 2] == "early"))
    expect_error(do.call(getVariogram, c(args, permutations = 9)))
    # split.by must not take the name of an output column
    tse$type <- tse$state
    args[[1]] <- tse
    args$split.by <- "type"
    expect_error(do.call(getVariogram, args), "output column")
})

test_that("Pairs are restricted to those from samples with subset.by TRUE", {
    data(hitchip1006, package = "miaTime")
    tse <- transformAssay(hitchip1006, method = "hellinger")
    tse$early <- tse$time < 3
    args <- list(tse, time.col = "time", group = "subject",
        assay.type = "hellinger", method = "euclidean", squared = TRUE,
        time.lag = c(1, 3))
    res <- do.call(getVariogram, c(args, subset.by = "early"))
    # The same as the local estimate on the pairs from early samples
    pairs <- miaTime:::.get_within_subject_pairs(tse, "time", "subject",
        "hellinger", "euclidean", NULL, TRUE, FALSE, list())
    first <- pmin(pairs$mid - pairs$h / 2)
    ref <- miaTime:::.local_variogram(pairs[first < 3, ], c(1, 3), 2, 0.95,
        c(NA, NA))
    expect_equal(res$estimate, ref$estimate)
    expect_equal(res$n.pairs, ref$n.pairs)
    expect_error(do.call(getVariogram, c(args, subset.by = "time")),
        "logical")
    tse$none <- FALSE
    args[[1]] <- tse
    expect_error(do.call(getVariogram, c(args, subset.by = "none")))
})

test_that("Pairs take the highest level of their two samples", {
    lev <- c("baseline", "exacerbation")
    a <- factor(c("baseline", "baseline", "exacerbation", NA, NA, NA),
        levels = lev)
    b <- factor(c("baseline", "exacerbation", "baseline", "exacerbation",
        "baseline", NA), levels = lev)
    expect_identical(as.character(miaTime:::.highest_level(a, b)),
        c("baseline", "exacerbation", "exacerbation", "exacerbation", NA, NA))
    # Every pair with an exacerbation sample, against pairs of two baseline
    # samples; samples of other states are kept for their pairs
    set.seed(1)
    n <- 160
    mat <- matrix(runif(5 * n), 5, dimnames = list(paste0("f", 1:5), NULL))
    tse <- SummarizedExperiment(assays = list(x = mat),
        colData = DataFrame(subject = rep(1:20, each = 8),
            time = rep(0:7, 20) + runif(n, 0, 0.5),
            state = factor(sample(c("baseline", "exacerbation", NA), n,
                replace = TRUE, prob = c(0.5, 0.3, 0.2)), levels = lev)))
    args <- list(tse, time.col = "time", group = "subject", assay.type = "x",
        method = "euclidean", squared = TRUE, time.lag = 2,
        split.by = "state")
    res <- do.call(getVariogram, c(args, pair.level = "highest"))
    pairs <- miaTime:::.get_within_subject_pairs(tse, "time", "subject", "x",
        "euclidean", "state", TRUE, FALSE, list(), "highest")
    st <- tse$state
    expect_true(all(pairs$level[st[pairs$k] %in% "exacerbation" |
        st[pairs$l] %in% "exacerbation"] == "exacerbation"))
    expect_true(all(pairs$level == "baseline" | st[pairs$k] %in%
        "exacerbation" | st[pairs$l] %in% "exacerbation"))
    expect_identical(unique(res$state), c("baseline", "exacerbation",
        "exacerbation / baseline"))
    expect_error(do.call(getVariogram, c(args, pair.level = "later")))
})
