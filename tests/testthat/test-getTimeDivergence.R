test_that("getTimeDivergence", {
  data("hitchip1006")
  sub_hitchip <- hitchip1006[,800:1151]
  expect_s4_class(sub_hitchip, "TreeSummarizedExperiment")
  hitchipTime_subset <- getTimeDivergence(sub_hitchip, field = "subject",
                                     time_interval = 1,
                                     time_field = "time")

  expect_s4_class(hitchipTime_subset, "SummarizedExperiment")
  expect_equal(dim(hitchipTime_subset), c(130,352))


})
