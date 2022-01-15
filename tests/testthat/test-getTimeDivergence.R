test_that("getTimeDivergence", {
  data("hitchip1006")
  se <- hitchip1006
  # Subset to speed up computing
  # Just pick 4 subjects with 1-5 time points
  se <- se[, colData(hitchip1006)$subject %in% c("900", "934", "843", "875")]
  se2 <- getTimeDivergence(se, group = "subject",
                                     time_interval = 1,
                                     time_field = "time")

  # Input and output classes should match
  # expect_equal(class(se), class(se2))

})
