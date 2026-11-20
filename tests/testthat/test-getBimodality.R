# Tests for .calculate_skewness (references calculated with moments::skewness)
test_that(".calculate_skewness works as expected", {
    # Symmetrical data
    expect_equal(.calculate_skewness(c(1, 2, 3)), 0)
    set.seed(12346)
    expect_equal(.calculate_skewness(rnorm(10000)), 0, tolerance = 0.05)
    # Positively skewed data
    expect_equal(
        round(.calculate_skewness(c(1, 1, 1, 10, 10, 10, 10, 10, 10, 10)), 2),
        -0.87)
    # Negatively skewed data
    expect_equal(
        round(.calculate_skewness(c(10, 10, 10, 10, 1, 1, 1)), 2), -0.29)
    # Constant data
    expect_equal(.calculate_skewness(rep(5, 10)), NaN)
})

# Tests for .calculate_kurtosis (references calculated with moments::kurtosis)
test_that(".calculate_kurtosis works as expected", {
    # Normalish data
    expect_equal(
        round(.calculate_kurtosis(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)), 2), 1.78)
    # High kurtosis data
    expect_equal(
        round(.calculate_kurtosis(c(1, 1, 1, 1, 50, 50, 50, 50, 50, 50)), 2),
        1.17)
    # Constant data should result in NaN
    expect_equal(.calculate_kurtosis(rep(5, 10)), NaN)
})

# Tests for .calculate_bimodality_coefficient (references calculated with
# mousetrap::bimodality_coefficient)
test_that(".calculate_bimodality_coefficient works as expected", {
    # Normalish data
    dat <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    skew <- .calculate_skewness(dat)
    kurt <- .calculate_kurtosis(dat)
    ref <- (1+skew^2)/(kurt+3)
    expect_equal(.calculate_bimodality_coefficient(dat), ref)
    # Constant data should result in NaN
    expect_equal(.calculate_bimodality_coefficient(rep(5, 10)), NaN)
})

# Create dummy data
tse <- makeTSE(nrow = 100, ncol = 20)
assayNames(tse) <- "counts"

# Tests for .calculate_bimodality
test_that(".calculate_bimodality works as expected", {
    result <- .calculate_bimodality(tse, assay.type = "counts", group = "group")
    expect_true(is.data.frame(result) || inherits(result, "DataFrame"))
    expect_equal(ncol(result), length(unique(tse$group)))
    expect_equal(nrow(result), nrow(tse))
})

# Tests for getBimodality
test_that("getBimodality works as expected", {
    result <- getBimodality(tse, assay.type = "counts")
    expect_true(is.data.frame(result) || inherits(result, "DataFrame"))
    expect_equal(ncol(result), 1)
    expect_equal(nrow(result), nrow(tse))
    ref <- .calculate_bimodality_coefficient(assay(tse, "counts")[1, ])
    expect_equal(result[[1]][[1]], ref)
})

# Tests for addBimodality
test_that("addBimodality works as expected", {
    tse <- addBimodality(tse, assay.type = "counts", name = "testing")
    expect_true("testing" %in% colnames(rowData(tse)))
    result <- getBimodality(tse, assay.type = "counts", name = "testing")
    expect_equal(rowData(tse)[, "testing", drop = FALSE], result)
})

# Tests that error messages work
test_that("addBimodality: errors", {
    addBimodality(tse, assay.type = "test") |> expect_error()
    addBimodality(tse, group = "test") |> expect_error()
    addBimodality(tse, name = 1) |> expect_error()
    addBimodality(tse, name = c("test", "test2")) |> expect_error()
    addBimodality(tse, name = TRUE) |> expect_error()
})

test_that("potential analysis gives the multimodality score of microbiome", {
    data(SilvermanAGutData)
    tse <- SilvermanAGutData[, colSums(assay(SilvermanAGutData, "counts")) > 0]
    tse <- transformAssay(tse, method = "relabundance")
    a <- assay(tse, "relabundance")
    # microbiome::multimodality(x, peak.threshold = 0, bs.iter = 30) with
    # set.seed(7) gives these scores
    expected <- c(seq_3 = 7 / 30, seq_5 = 9 / 30, seq_9 = 18 / 30)
    for( feature in names(expected) ){
        set.seed(7)
        score <- .calculate_multimodality(as.numeric(a[feature, ]),
            peak.threshold = 0, bs.iter = 30)
        expect_equal(score, expected[[feature]])
    }
    # Modes where shallow maxima are dropped and neighbours merged; the
    # positions are those of microbiome:::potential_univariate()
    expect_equal(.potential_maxima(as.numeric(a["seq_16", ]), 1, 1, 1),
        c(4, 30, 65, 99))
    expect_equal(.potential_maxima(as.numeric(a["seq_1", ]), 1, 1, 1),
        c(12, 40))
})

test_that("getBimodality(method = 'potential') separates one and two modes", {
    set.seed(3)
    x <- rbind(
        one = rnorm(200),
        two = c(rnorm(100, mean = 0), rnorm(100, mean = 6)),
        flat = rep(1, 200))
    se <- SummarizedExperiment(assays = list(counts = x),
        colData = DataFrame(g = rep(c("a", "b"), 100)))
    res <- getBimodality(se, method = "potential", bs.iter = 50)
    expect_s4_class(res, "DataFrame")
    expect_identical(rownames(res), rownames(x))
    expect_lt(res["one", 1], 0.3)
    expect_gt(res["two", 1], 0.9)
    expect_equal(res["flat", 1], 0)
    # Missing values are ignored
    x_na <- x
    x_na["two", 1:5] <- NA
    se_na <- SummarizedExperiment(assays = list(counts = x_na))
    expect_false(anyNA(getBimodality(se_na, method = "potential",
        bs.iter = 5)[[1]]))
    # One column per group, and stored in rowData by addBimodality
    res_g <- getBimodality(se, method = "potential", bs.iter = 10, group = "g")
    expect_identical(colnames(res_g), c("bimodality_a", "bimodality_b"))
    se <- addBimodality(se, method = "potential", bs.iter = 10, name = "pot")
    expect_true("pot" %in% colnames(rowData(se)))
    expect_error(getBimodality(se, method = "other"))
})

test_that("a higher peak.threshold keeps fewer modes", {
    set.seed(5)
    x <- c(rnorm(100), rnorm(100, mean = 3))
    n_modes <- vapply(c(0, 1, 10), function(t){
        length(.potential_maxima(x, t, 1, 1))
    }, numeric(1))
    expect_identical(n_modes, c(2, 1, 1))
})
