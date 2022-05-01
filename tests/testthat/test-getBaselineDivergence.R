test_that("getBaselineDivergence", {
  data("hitchip1006")
  tse <- hitchip1006
  # Subset to speed up computing
  # Just pick 4 subjects with 1-5 time points
  tse <- tse[, colData(tse)$subject %in% c("900", "934", "843", "875", "836")]
  tse2 <- getBaselineDivergence(tse, group = "subject", time_field = "time")

  # Input and output classes should match
  expect_equal(class(tse), class(tse2))

  # A subject to check time difference calculation
  time2 <- colData(tse2)[, "time"][which(colData(tse2)[, "subject"] == "843")]
  time_diff_2 <- colData(tse2)[, "time_from_baseline"][which(colData(tse2)[, "subject"] == "843")]
  expect_true(all(time2==time_diff_2))

  # Test divergences
  inds0 <- which(colData(tse)[, "subject"] == "843")  
  inds <- which(colData(tse2)[, "subject"] == "843")
  original.divergence <- as.matrix(vegan::vegdist(t(assay(tse[, inds0], "counts"))))[,1]
  calculated.divergence <- colData(tse2)[inds, "divergence_from_baseline"]
  expect_true(all(original.divergence==calculated.divergence))

  # Should also work when baseline is not 0  
  inds <- which(colData(tse)[, "subject"] == "843")[2:5]
  tse2 <- getBaselineDivergence(tse[, inds], group = "subject", time_field = "time")
  time2 <- colData(tse[, inds])[, "time"] - min(colData(tse[, inds])[, "time"])
  time_diff_2 <- colData(tse2)[, "time_from_baseline"]
  expect_true(all(time2==time_diff_2))


})
