test_that("getStepwiseDivergence", {
  data(hitchip1006)
  tse <- hitchip1006
  # Subset to speed up computing
  # Just pick 4 subjects with 1-5 time points
  tse <- tse[, colData(tse)$subject %in% c("900", "934", "843", "875")]
  tse2 <- addStepwiseDivergence(tse, group = "subject",
                                time_interval = 1,
                                time_field = "time",
                                assay.type="counts",
                                FUN = vegan::vegdist,
                                method="bray")

  # Trying to add new coldata field with the same name
  expect_warning(tse2 <- addStepwiseDivergence(tse2, group = "subject",
                                     time_interval = 1,
                                     time_field = "time"))

  # Input and output classes should match
  expect_equal(class(tse), class(tse2))

  # A subject to check time difference calculation
  obs_diff <- colData(tse2)[
      which(colData(tse2)[, "subject"] == "843"), "time_difference"]
  exp_diff <- c(NA,diff(colData(tse)[
      which(colData(tse)[, "subject"] == "843"), "time"]))
  expect_equal(obs_diff, exp_diff)

  # n > 1
  tse3 <- addStepwiseDivergence(tse, group = "subject",
                                time_interval = 2,
                                time_field = "time",
                                assay.type="counts",
                                FUN = vegan::vegdist,
                                method="bray")
  
  time_invertal <- 2

  time3 <- colData(tse3)[, "time"][which(colData(tse3)[, "subject"] == "843")]

  time_dif_3 <- colData(tse3)[, "time_difference"][
      which(colData(tse3)[, "subject"] == "843")]

  # number of divergences (n-k) check
  divergence_number <- length(time3) - time_invertal

  divergence_calculated <- length(which(!is.na(time_dif_3) == TRUE))

  expect_equal(divergence_number, divergence_calculated)

  # interval check
  calculated_diff <- time3[(1+ 2):length(time3)] - 
      time3[seq_len(length(time3)-2)]

  manual_diff <- c(rep(NA, length(time3) - 
      length(calculated_diff)), calculated_diff)

  expect_equal(time_dif_3, manual_diff)

  # object with single time point has NA instead of divergence values
  sub_hitchip <- hitchip1006[, 
      colData(hitchip1006)$subject %in% c("900","843", "139")]
  subset <- addStepwiseDivergence(sub_hitchip, group = "subject",
                                time_interval = 1,
                                time_field = "time",
                                assay.type="counts",
                                FUN = vegan::vegdist,
                                method="bray")

  expect_true(all(is.na(colData(subset)[, "time_divergence"][
      which(duplicated(colData(subset)[, "subject"]) == FALSE)])))


  # Test vegan distances
  tse2 <- addStepwiseDivergence(tse, group = "subject",
                                  time_interval = 1,
                                  time_field = "time",
                                  assay.type="counts",
                                  FUN = vegan::vegdist,
                                  method="bray",
                                  name_timedifference="timedifference",
                                  name_divergence="timedivergence")
  
  # Test vegan distances
  tse2 <- addStepwiseDivergence(tse2, group = "subject",
                                time_interval = 1,
                                time_field = "time",
                                assay.type="counts",
                                FUN = vegan::vegdist,
                                method="euclidean",
                                name_timedifference="timedifference2",
                                name_divergence="timedivergence2")
				     
  # Time differences should still match
  expect_true(identical(tse2$timedifference, tse2$timedifference2))
  # ... but divergences should be different (bray vs. euclid)
  expect_true(!identical(tse2$timedivergence, tse2$timedivergence2))
  
  ## Test with ordination values
  tse2 <- scater::runMDS(tse2, FUN = vegan::vegdist, method = "bray",
                        name = "PCoA_BC", exprs_values = "counts",
                        na.rm = TRUE, ncomponents=4)
  # testing with all ordination components; n_dimred=NULL --> all 4 components
  tse2 <- addStepwiseDivergence(tse2, group = "subject",
                                time_interval = 1,
                                time_field = "time",
                                name_timedifference="timedifference_ord_4",
                                name_divergence="timedivergence_ord_4",
                                dimred = "PCoA_BC",
                                FUN=vegan::vegdist,
                                method="euclidean")
  # Time differences should still match
  expect_true(identical(tse2$timedifference_ord_4, tse2$timedifference))
  # ordination based divergence values should not be equal to the ones on counts
  expect_true(!identical(tse2$timedivergence_ord_4, tse2$timedivergence))
  
  # testing with 2 ordination components
  tse2 <- addStepwiseDivergence(tse2, group = "subject",
                                time_interval = 1,
                                time_field = "time",
                                name_timedifference="timedifference_ord_2",
                                name_divergence="timedivergence_ord_2",
                                dimred = "PCoA_BC",
                                n_dimred = 2,
                                FUN=vegan::vegdist,
                                method="euclidean")
  # Time differences should still match
  expect_true(identical(tse2$timedifference_ord_2, tse2$timedifference_ord_4))
  # not same values as using 4 components
  expect_true(!identical(tse2$timedivergence_ord_2, tse2$timedivergence_ord_4))
  
  ## testing with altExp
  SingleCellExperiment::altExp(tse2, "Family") <- mia::agglomerateByRank(tse2, 
                                                        rank="Family")
  tse2 <- addStepwiseDivergence(tse2, group = "subject",
                                time_interval = 1,
                                time_field = "time",
                                altexp="Family",
                                name_timedifference="timedifference_Fam",
                                name_divergence="timedivergence_Fam")
  # Time differences should still match
  expect_true(identical(tse2$timedifference_Fam, tse2$timedifference))
  # divergence values based on Family rank counts should not be equal to the
  # ones with Genus counts
  expect_true(!identical(tse2$timedivergence_Fam, tse2$timedivergence))
})
