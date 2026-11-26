data(minimalgut)
tse <- minimalgut

test_that("getRankShift checks its input", {
    expect_error(getRankShift(tse, time.col = "missing"))
    expect_error(getRankShift(tse, time.col = "StudyIdentifier"), "numeric")
    expect_error(getRankShift(tse, time.col = "Time_hr_num",
        assay.type = "missing"))
    expect_error(addRankShift(tse, time.col = "Time_hr_num", name = "a"))
})

test_that("getRankShift equals codyn::rank_shift", {
    res <- getRankShift(tse, time.col = "Time_hr_num",
        group = "StudyIdentifier")
    a <- tse$StudyIdentifier == "Bioreactor A"
    shift <- res[["rank_shift"]][a][order(tse$Time_hr_num[a])]
    # codyn::rank_shift() of Bioreactor A, without zero abundances
    expect_equal(shift[2:5], c(46, 22, 24, 18) / 15)
    expect_true(is.na(shift[1L]))
    expect_identical(colnames(res), c("rank_shift", "time_diff", "ref_samples"))
    expect_identical(rownames(res), colnames(tse))
})

test_that("getRankShift compares each sample with the previous one", {
    res <- getRankShift(tse, time.col = "Time_hr_num",
        group = "StudyIdentifier")
    expect_equal(sum(is.na(res[["rank_shift"]])), 3L)
    t0 <- tse$Time_hr_num == 0
    expect_true(all(is.na(res[["rank_shift"]][t0])))
    later <- which(!t0)[1L]
    ref <- res[["ref_samples"]][later]
    expect_equal(res[["time_diff"]][later],
        tse$Time_hr_num[later] - tse[, ref]$Time_hr_num)
    # time.lag picks the sample closest to the given lag instead
    res_lag <- getRankShift(tse, time.col = "Time_hr_num",
        group = "StudyIdentifier", time.lag = 24)
    expect_false(identical(res_lag[["ref_samples"]], res[["ref_samples"]]))
})

test_that(".rank_shift is the mean absolute change in rank", {
    expect_equal(.rank_shift(c(1, 2, 3), c(10, 20, 30)), 0)
    expect_equal(.rank_shift(c(1, 2, 3), c(3, 2, 1)), 4 / 3)
    # Features absent from either sample are left out
    expect_equal(.rank_shift(c(0, 1, 2), c(5, 2, 1)), 1)
    expect_true(is.na(.rank_shift(c(0, 1), c(1, 0))))
})

test_that("addRankShift adds the results to colData", {
    x <- addRankShift(tse, time.col = "Time_hr_num", group = "StudyIdentifier",
        name = c("shift", "dt", "prev"))
    expect_true(all(c("shift", "dt", "prev") %in% colnames(colData(x))))
    expect_equal(x[["shift"]], getRankShift(tse, time.col = "Time_hr_num",
        group = "StudyIdentifier")[["rank_shift"]])
})
