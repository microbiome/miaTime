test_that("getVariogramContribution input checks", {
    data(hitchip1006, package = "miaTime")
    tse <- hitchip1006
    expect_error(getVariogramContribution(tse, time.col = "x",
        group = "subject"))
    expect_error(getVariogramContribution(tse, time.col = "time",
        group = "subject", assay.type = "x"))
    expect_error(getVariogramContribution(tse, time.col = "time",
        group = "subject", time.lag = -1))
    expect_error(getVariogramContribution(tse, time.col = "time",
        group = "subject", time.lag = 1, window = 1))
    expect_error(getVariogramContribution(tse, time.col = "time",
        group = "subject", time.lag = 1, conf.level = 2))
})

test_that("Contributions split D(L) exactly", {
    data(hitchip1006, package = "miaTime")
    tse <- transformAssay(hitchip1006, method = "hellinger")
    args <- list(tse, time.col = "time", group = "subject",
        assay.type = "hellinger", time.lag = c(1, 3))
    res <- do.call(getVariogramContribution, args)
    expect_s4_class(res, "DataFrame")
    expect_equal(nrow(res), 2 * nrow(tse))
    expect_setequal(res$feature, rownames(tse))
    D <- do.call(getVariogram, c(args, method = "euclidean", squared = TRUE))
    sums <- tapply(res$estimate, res$lag, sum)
    expect_equal(as.numeric(sums), D$estimate)
    expect_equal(unname(tapply(res$share, res$lag, sum))[D$supported], 1)
    # With one feature, the contribution is D(L) with the same standard error
    f <- res$feature[which.max(res$estimate)]
    one <- do.call(getVariogram, c(list(tse[f, ]), args[-1],
        method = "euclidean", squared = TRUE))
    r <- res[res$feature == f, ]
    expect_equal(r$estimate, one$estimate)
    expect_equal(r$se, one$se * one$estimate)
    # At a stated time
    res <- do.call(getVariogramContribution, c(args, time.point = 3))
    D <- do.call(getVariogram, c(args, method = "euclidean", squared = TRUE,
        time.point = 3))
    expect_equal(as.numeric(tapply(res$estimate, res$lag, sum)), D$estimate)
    expect_true(all(res$time.point == 3))
})

test_that("Contributions are compared between levels", {
    data(hitchip1006, package = "miaTime")
    tse <- transformAssay(hitchip1006, method = "hellinger")
    args <- list(tse, time.col = "time", group = "subject",
        assay.type = "hellinger", time.lag = 2.4, split.by = "sex")
    res <- do.call(getVariogramContribution, args)
    expect_identical(unique(res$sex), c("female", "male", "male - female"))
    diff <- res[res$type == "difference", ]
    expect_equal(diff$estimate,
        res$estimate[res$sex == "male"] - res$estimate[res$sex == "female"])
    expect_true(all(is.na(diff$share)))
    expect_true(all(diff$supported))
    # Features that never change have no p-value
    expect_true(all(diff$p.value >= 0 & diff$p.value <= 1, na.rm = TRUE))
    expect_identical(is.na(diff$p.value), diff$se == 0)
    # Sample-level state: the differences sum to b1 of the regression on d.
    # Twenty subjects with eight samples each, alternating between states
    set.seed(1)
    n <- 160
    mat <- matrix(runif(5 * n), 5, dimnames = list(paste0("f", 1:5), NULL))
    tse <- SummarizedExperiment(assays = list(x = mat),
        colData = DataFrame(subject = rep(1:20, each = 8),
            time = rep(0:7, 20) + runif(n, 0, 0.5),
            state = rep(c("a", "b"), 80)))
    res <- suppressMessages(getVariogramContribution(tse, time.col = "time",
        group = "subject", assay.type = "x", time.lag = 2,
        split.by = "state"))
    pairs <- suppressMessages(miaTime:::.get_within_subject_pairs(tse,
        "time", "subject", "x", "euclidean", "state", TRUE, FALSE, list()))
    f <- miaTime:::.pair_level_design(pairs, "b", 2, NA_real_, 2)
    b <- f$Ai %*% crossprod(f$X, f$w * pairs$d[f$rows])
    diff <- res[res$type == "difference", ]
    expect_true(all(diff$supported))
    expect_equal(sum(diff$estimate), b[[2]])
})
