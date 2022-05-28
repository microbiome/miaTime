test_that("getBaselineDivergence", {

  library(dplyr)
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

  # -----------------------------------------------------------

  # devtools::load_all("~/Rpackages/microbiome/miaverse/miaTime/")

  data("hitchip1006")
  tse <- hitchip1006
  # Just pick 1 subject with many time points
  tse <- tse[, colData(tse)$subject == "843"] # The baseline time point 0 is Sample-843
  
  # Should now work also without the "group" argument because there is just a single group
  tse2a <- getBaselineDivergence(tse, time_field = "time")
  tse2b <- getBaselineDivergence(tse, group=1, time_field = "time")
  expect_identical(tse2a, tse2b)

  # Define the baseline sample manually
  tse2c <- getBaselineDivergence(tse, time_field = "time", baseline_sample="Sample-843")
  tse2d <- getBaselineDivergence(tse, time_field = "time", baseline_sample="Sample-1075")
  # Now the times from baseline should be shifted and dissimilarities differ

  # Sample baseline when the zero time baseline is automatically checked or manually set
  expect_true(all(tse2b$time_from_baseline==tse2c$time_from_baseline))
  # The shifted case (different, middle sample as baseline)
  expect_true(all(tse2c$time_from_baseline == tse2d$time_from_baseline + 0.7))

  tse <- hitchip1006
  # Subset to speed up computing
  # Just pick 4 subjects with 1-5 time points
  tse <- tse[, colData(tse)$subject %in% c("900", "934", "843", "875", "836")]
  tse2e <- getBaselineDivergence(tse[, colData(tse)$subject == "843"], group="subject", time_field = "time")      
  tse2f <- getBaselineDivergence(tse, group = "subject", time_field = "time")
  tse2g <- getBaselineDivergence(tse, group = "subject", time_field = "time", baseline_sample="Sample-1075")  
  expect_identical(colData(tse2e)["Sample-843", "time_from_baseline"], colData(tse2f)["Sample-843", "time_from_baseline"])
  expect_identical(colData(tse2e)["Sample-843", "time_from_baseline"] - 0.7, colData(tse2g)["Sample-843", "time_from_baseline"])  

  # Test with full baseline list
  baselines <- c("Sample-1041", "Sample-1075",  "Sample-875", "Sample-900", "Sample-934")
  names(baselines) <- names(split(colnames(tse), as.character(tse$subject)))
  tse2h <- getBaselineDivergence(tse, group = "subject", time_field = "time", baseline_sample=baselines)
  expect_identical(colData(tse2h)["Sample-843", "time_from_baseline"], colData(tse2g)["Sample-843", "time_from_baseline"])    

  # Single baseline
  tse2i <- getBaselineDivergence(tse, group = "subject", time_field = "time", baseline_sample=tse[, "Sample-1075"])
  expect_identical(colData(tse2i)["Sample-1075", "time_from_baseline"], colData(tse2g)["Sample-1075", "time_from_baseline"])
  expect_identical(colData(tse2i)["Sample-843", "time_from_baseline"] + 0.7, colData(tse2g)["Sample-1075", "time_from_baseline"])  
  
})
