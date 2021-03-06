test_that("getStepwiseDivergence", {
  data("hitchip1006")
  tse <- hitchip1006
  # Subset to speed up computing
  # Just pick 4 subjects with 1-5 time points
  tse <- tse[, colData(tse)$subject %in% c("900", "934", "843", "875")]
  tse2 <- getStepwiseDivergence(tse, group = "subject",
                                     time_interval = 1,
                                     time_field = "time")

  # Input and output classes should match
  expect_equal(class(tse), class(tse2))

  # A subject to check time difference calculation
  time2 <- colData(tse2)[, "time"][which(colData(tse2)[, "subject"] == "843")]

  time_diff_2 <- colData(tse2)[, "time_difference"][which(colData(tse2)[, "subject"] == "843")]

  diff_vector <- c(NA,diff(time2))

  expect_equal(time_diff_2, diff_vector)

  # n > 1
  tse3 <- getStepwiseDivergence(tse, group = "subject",
                           time_interval = 2,
                           time_field = "time")
  time_invertal <- 2

  time3 <- colData(tse3)[, "time"][which(colData(tse3)[, "subject"] == "843")]

  time_dif_3 <- colData(tse3)[, "time_difference"][which(colData(tse3)[, "subject"] == "843")]

  # number of divergences (n-k) check
  divergence_number <- length(time3) - time_invertal

  divergence_calculated <- length(which(!is.na(time_dif_3) == TRUE))

  expect_equal(divergence_number, divergence_calculated)

  # interval check
  calculated_diff <- time3[(1+ 2):length(time3)] - time3[seq_len(length(time3)-2)]

  manual_diff <- c(rep(NA, length(time3) - length(calculated_diff)), calculated_diff)

  expect_equal(time_dif_3, manual_diff)

  # object with single time point has NA instead of divergence values
  sub_hitchip <- hitchip1006[, colData(hitchip1006)$subject %in% c("900","843", "139")]
  subset <- getStepwiseDivergence(sub_hitchip, group = "subject",
                                      time_interval = 1,
                                      time_field = "time")

  expect_true(all(is.na(colData(subset)[, "time_divergence"][which(duplicated(colData(subset)[, "subject"]) == FALSE)])))

})
