test_that("getBaselineDivergence", {

  library(dplyr)
  data(hitchip1006)
  tse <- hitchip1006
  # Subset to speed up computing
  # Just pick 4 subjects with 1-5 time points
  tse <- tse[, colData(tse)$subject %in% c("900", "934", "843", "875", "836")]
  tse2 <- addBaselineDivergence(tse, group = "subject", time.col = "time", 
              name.time = "time_from_baseline", 
              name = "divergence_from_baseline")

  # Input and output classes should match
  expect_equal(class(tse), class(tse2))

  # A subject to check time difference calculation
  time2 <- colData(tse2)[, "time"][which(colData(tse2)[, "subject"] == "843")]
  time_diff_2 <- colData(tse2)[, "time_from_baseline"][
    which(colData(tse2)[, "subject"] == "843")]
  expect_true(all(time2==time_diff_2))

  # Test divergences
  inds0 <- which(colData(tse)[, "subject"] == "843")  
  inds <- which(colData(tse2)[, "subject"] == "843")
  original.divergence <- as.matrix(
    vegan::vegdist(t(assay(tse[, inds0], "counts"))))[,1]
  calculated.divergence <- colData(tse2)[inds, "divergence_from_baseline"]
  expect_true(all(original.divergence==calculated.divergence))

  # Should also work when baseline is not 0  
  inds <- which(colData(tse)[, "subject"] == "843")[2:5]
  tse2 <- addBaselineDivergence(tse[, inds], group = "subject", 
              time.col = "time",
              name.time = "time_from_baseline", 
              name = "divergence_from_baseline")
  time2 <- colData(tse[, inds])[, "time"] - min(colData(tse[, inds])[, "time"])
  time_diff_2 <- colData(tse2)[, "time_from_baseline"]
  expect_true(all(time2==time_diff_2))

  # -----------------------------------------------------------

  # devtools::load_all("~/Rpackages/microbiome/miaverse/miaTime/")

  data(hitchip1006)
  tse <- hitchip1006
  # Just pick 1 subject with many time points
  # The baseline time point 0 is Sample-843
  tse <- tse[, colData(tse)$subject == "843"] 
  
  tse2b <- addBaselineDivergence(tse, group="subject", time.col = "time")
  # Define the baseline sample manually
  tse2c <- addBaselineDivergence(tse, time.col = "time", group="subject",
               baseline_sample="Sample-843",
               name.time = "time_from_baseline", 
               name = "divergence_from_baseline")
  tse2d <- addBaselineDivergence(tse, time.col = "time", group="subject",
               baseline_sample="Sample-1075",
               name.time = "time_from_baseline", 
               name = "divergence_from_baseline")
  # Now the times from baseline should be shifted and dissimilarities differ

  # Sample baseline when the zero time baseline is automatically checked or 
  # manually set
  expect_true(all(tse2b$time_from_baseline==tse2c$time_from_baseline))
  # The shifted case (different, middle sample as baseline)
  expect_true(all(tse2c$time_from_baseline == tse2d$time_from_baseline + 0.7))

  tse <- hitchip1006
  # Subset to speed up computing
  # Just pick 4 subjects with 1-5 time points
  tse <- tse[, colData(tse)$subject %in% c("900", "934", "843", "875", "836")]
  tse2e <- addBaselineDivergence(tse[, colData(tse)$subject == "843"], 
               group="subject", time.col = "time",
               name.time = "time_from_baseline", 
               name = "divergence_from_baseline")      
  tse2f <- addBaselineDivergence(tse, group = "subject", time.col = "time", 
               name.time = "time_from_baseline", 
               name = "divergence_from_baseline")
  tse2g <- addBaselineDivergence(tse, group = "subject", time.col = "time", 
               baseline_sample="Sample-1075", 
               name.time = "time_from_baseline", 
               name = "divergence_from_baseline")  
  expect_identical(colData(tse2e)["Sample-843", "time_from_baseline"], 
                   colData(tse2f)["Sample-843", "time_from_baseline"])
  expect_identical(colData(tse2e)["Sample-843", "time_from_baseline"] - 0.7, 
                   colData(tse2g)["Sample-843", "time_from_baseline"])  

  # Test with full baseline list
  baselines <- c("Sample-1041", "Sample-1075",  
                 "Sample-875", "Sample-900", "Sample-934")
  names(baselines) <- names(split(colnames(tse), as.character(tse$subject)))
  tse2h <- addBaselineDivergence(tse, group = "subject", time.col = "time", 
               baseline_sample=baselines, 
               name.time = "time_from_baseline", 
               name = "divergence_from_baseline")
  expect_identical(colData(tse2h)["Sample-843", "time_from_baseline"], 
                   colData(tse2g)["Sample-843", "time_from_baseline"])    

  # Single baseline
  tse2i <- addBaselineDivergence(tse, group = "subject", time.col = "time", 
               baseline_sample=tse[, "Sample-1075"], 
               name.time = "time_from_baseline", 
               name = "divergence_from_baseline")
  expect_identical(colData(tse2i)["Sample-1075", "time_from_baseline"], 
                   colData(tse2g)["Sample-1075", "time_from_baseline"])
  expect_identical(colData(tse2i)["Sample-843", "time_from_baseline"] + 0.7, 
                   colData(tse2g)["Sample-1075", "time_from_baseline"])  
  
  ## Test with ordination values
  tse <- scater::runMDS(tse, FUN = vegan::vegdist, method = "bray",
                         name = "PCoA_BC", exprs_values = "counts",
                         na.rm = TRUE, ncomponents=4)
  # testing with all ordination components; ndimred=NULL --> all 4 components
  tse2 <- addBaselineDivergence(tse, group = "subject",
                                time.col = "time",
                                name.time="time_from_baseline_ord_4",
                                name="divergence_from_baseline_ord_4",
                                dimred = "PCoA_BC",
                                dis.fun=vegan::vegdist,
                                method="euclidean")
  # Time differences should still match
  expect_true(identical(tse2$time_from_baseline_ord_4, 
                        tse2f$time_from_baseline))
  # ordination based divergence values should not be equal to the ones on counts
  expect_false(identical(tse2$divergence_from_baseline_ord_4, 
                       tse2f$divergence_from_baseline))
  # testing with 2 ordination components
  tse2 <- addBaselineDivergence(tse2, group = "subject",
                                time.col = "time",
                                name.time="time_from_baseline_ord_2",
                                name="divergence_from_baseline_ord_2",
                                dimred = "PCoA_BC",
                                ndimred = 2,
                                dis.fun=vegan::vegdist,
                                method="euclidean")
  # Time differences should still match
  expect_true(identical(tse2$time_from_baseline_ord_4, 
                        tse2$time_from_baseline_ord_2))
  # ordination based divergence values should not be equal to the ones on counts
  expect_false(identical(tse2$divergence_from_baseline_ord_4, 
                         tse2$divergence_from_baseline_ord_2))
  ## testing with altExp
  SingleCellExperiment::altExp(tse2, "Family") <- agglomerateByRank(tse2, 
                                                      rank="Family")
  tse2 <- addBaselineDivergence(tse2, group = "subject",
                                time.col = "time",
                                altexp="Family",
                                name.time="time_from_baseline_Fam",
                                name="divergence_from_baseline_Fam")
  # Time differences should still match
  expect_true(identical(tse2$time_from_baseline_Fam, tse2f$time_from_baseline))
  # divergence values based on Family rank counts should not be equal to the
  # ones with Genus counts
  expect_false(identical(tse2$divergence_from_baseline_Fam, 
                         tse2f$divergence_from_baseline))
})
