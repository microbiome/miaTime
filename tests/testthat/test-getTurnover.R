data(HMP_2019_ibdmdb)
tse <- HMP_2019_ibdmdb
time <- "days_from_first_collection"

test_that("getTurnover checks its input", {
    expect_error(getTurnover(tse, time.col = "missing"))
    expect_error(getTurnover(tse, time.col = "subject_id"), "numeric")
    expect_error(getTurnover(tse, time.col = time, assay.type = "missing"))
    expect_error(getTurnover(tse, time.col = time,
        assay.type = "relabundance", method = "missing"))
    expect_error(getTurnover(tse, time.col = time,
        assay.type = "relabundance", detection = "a"))
    expect_error(addTurnover(tse, time.col = time,
        assay.type = "relabundance", name = "a"))
})

test_that("getTurnover equals codyn::turnover", {
    # codyn::turnover() of subjects C3003 and C3005, without zero abundances
    get <- function(method){
        res <- getTurnover(tse, time.col = time, group = "subject_id",
            assay.type = "relabundance", method = method)
        res <- res[!is.na(res[["turnover"]]), ]
        res[["turnover"]][match(c("C3003", "C3005"),
            tse[, rownames(res)]$subject_id)]
    }
    expect_equal(get("total"), c(23 / 67, 46 / 87))
    expect_equal(get("appearance"), c(12 / 67, 9 / 87))
    expect_equal(get("disappearance"), c(11 / 67, 37 / 87))
})

test_that("getTurnover compares each sample with the previous one", {
    res <- getTurnover(tse, time.col = time, group = "subject_id",
        assay.type = "relabundance")
    expect_identical(colnames(res), c("turnover", "time_diff", "ref_samples"))
    expect_identical(rownames(res), colnames(tse))
    # Same reference samples as stepwise divergence
    div <- getStepwiseDivergence(tse, time.col = time, group = "subject_id",
        assay.type = "relabundance")
    expect_identical(res[["ref_samples"]], div[["ref_samples"]])
    expect_identical(is.na(res[["turnover"]]), is.na(div[["divergence"]]))
    expect_true(all(res[["turnover"]] >= 0 & res[["turnover"]] <= 1,
        na.rm = TRUE))
})

test_that("a higher detection threshold makes more features absent", {
    res <- vapply(c(0, 0.01, 0.1), function(d)
        getTurnover(tse, time.col = time, group = "subject_id",
            assay.type = "relabundance", detection = d)[["turnover"]] |>
            mean(na.rm = TRUE), numeric(1))
    expect_false(any(duplicated(res)))
})

test_that(".turnover counts gains and losses among features present", {
    cur <- c(TRUE, TRUE, FALSE, FALSE)
    prev <- c(TRUE, FALSE, TRUE, FALSE)
    expect_equal(.turnover(cur, prev, "total"), 2 / 3)
    expect_equal(.turnover(cur, prev, "appearance"), 1 / 3)
    expect_equal(.turnover(cur, prev, "disappearance"), 1 / 3)
    expect_equal(.turnover(cur, cur, "total"), 0)
    expect_true(is.na(.turnover(c(FALSE, FALSE), c(FALSE, FALSE), "total")))
})

test_that("addTurnover adds the results to colData", {
    x <- addTurnover(tse, time.col = time, group = "subject_id",
        assay.type = "relabundance", method = "appearance",
        name = c("gained", "dt", "prev"))
    expect_true(all(c("gained", "dt", "prev") %in% colnames(colData(x))))
    expect_equal(x[["gained"]], getTurnover(tse, time.col = time,
        group = "subject_id", assay.type = "relabundance",
        method = "appearance")[["turnover"]])
})
