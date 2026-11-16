data(minimalgut)
tse <- transformAssay(minimalgut, method = "relabundance")
m <- assay(tse, "relabundance")

test_that("getFeatureIndices checks its input", {
    expect_error(getFeatureIndices(tse, assay.type = "missing"))
    expect_error(getFeatureIndices(tse, assay.type = "relabundance",
        index = "missing"))
    expect_error(getFeatureIndices(tse, assay.type = "relabundance",
        index = c("cv", "cv")), "duplicates")
    expect_error(getFeatureIndices(tse, assay.type = "relabundance",
        group = "missing"))
    expect_error(addFeatureIndices(tse, assay.type = "relabundance",
        index = c("cv", "variance"), name = "a"))
    # Stability needs the time column
    expect_error(getFeatureIndices(tse, assay.type = "relabundance",
        index = "stability"))
})

test_that("getFeatureIndices equals direct calculations", {
    res <- getFeatureIndices(tse, assay.type = "relabundance")
    expect_identical(colnames(res),
        c("variance", "cv", "skewness", "kurtosis", "bimodality"))
    expect_identical(rownames(res), rownames(tse))
    expect_equal(res[["variance"]], unname(apply(m, 1, var)))
    expect_equal(res[["cv"]], unname(apply(m, 1, sd) / rowMeans(m)))
    z <- (m - rowMeans(m)) / sqrt(rowMeans((m - rowMeans(m))^2))
    expect_equal(res[["skewness"]], unname(rowMeans(z^3)))
    expect_equal(res[["kurtosis"]], unname(rowMeans(z^4)))
    expect_equal(res[["bimodality"]], getBimodality(tse,
        assay.type = "relabundance")[["bimodality"]])
})

test_that("getFeatureIndices works within groups", {
    res <- getFeatureIndices(tse, assay.type = "relabundance",
        index = c("variance", "stability"), group = "StudyIdentifier",
        time.col = "Time.hr")
    groups <- unique(tse$StudyIdentifier)
    expect_setequal(colnames(res), c(paste0("variance_", groups),
        paste0("stability_", groups)))
    a <- tse$StudyIdentifier == "Bioreactor A"
    expect_equal(res[["variance_Bioreactor A"]],
        unname(apply(m[, a], 1, var)))
    stab <- getStability(tse, assay.type = "relabundance",
        group = "StudyIdentifier", time.col = "Time.hr")
    expect_equal(res[["stability_Bioreactor A"]],
        stab[["stability_Bioreactor A"]])
})

test_that("addFeatureIndices adds the indices to rowData", {
    x <- addFeatureIndices(tse, assay.type = "relabundance",
        index = c("cv", "bimodality"), name = c("cv_all", "bc"))
    expect_true(all(c("cv_all", "bc") %in% colnames(rowData(x))))
    expect_equal(rowData(x)[["cv_all"]], getFeatureIndices(tse,
        assay.type = "relabundance", index = "cv")[["cv"]])
})
