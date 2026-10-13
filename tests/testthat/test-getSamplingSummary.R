test_that("getSamplingSummary summarizes each group", {
    tse <- SummarizedExperiment(
        assays = list(counts = matrix(1, nrow = 2L, ncol = 6L)),
        colData = DataFrame(
            time = c(0, 1, 1, 4, 2, 5),
            subject = c("a", "a", "a", "a", "b", "c")))
    tse[["time"]][[6L]] <- NA

    res <- getSamplingSummary(tse, time.col = "time", group = "subject")

    expect_s4_class(res, "DataFrame")
    expect_equal(rownames(res), c("a", "b", "c"))
    expect_equal(res$n_samples, c(4L, 1L, 1L))
    expect_equal(res$n_time_points, c(3L, 1L, 0L))
    expect_equal(res$start, c(0, 2, NA))
    expect_equal(res$end, c(4, 2, NA))
    expect_equal(res$span, c(4, 0, NA))
    expect_equal(res$median_gap, c(2, NA, NA))
    expect_equal(res$max_gap, c(3, NA, NA))
})

test_that("getSamplingSummary without group summarizes all samples", {
    data(hitchip1006)
    res <- getSamplingSummary(hitchip1006, time.col = "time")

    expect_equal(nrow(res), 1L)
    expect_equal(res$n_samples, ncol(hitchip1006))
})

test_that("getSamplingSummary matches the subjects in the data", {
    data(hitchip1006)
    res <- getSamplingSummary(
        hitchip1006, time.col = "time", group = "subject")

    expect_equal(nrow(res), length(unique(hitchip1006$subject)))
    expect_equal(sum(res$n_samples), ncol(hitchip1006))
    expect_true(all(is.na(res$max_gap[res$n_time_points == 1L])))
})

test_that("getSamplingSummary checks input", {
    data(hitchip1006)
    expect_error(getSamplingSummary(hitchip1006, time.col = "subject"))
    expect_error(getSamplingSummary(hitchip1006, time.col = "missing"))
    expect_error(getSamplingSummary(
        hitchip1006, time.col = "time", group = "missing"))
})
