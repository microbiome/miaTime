

tse2 <- addBaselineDivergence(
    tse,
    group = "subject",
    time_field = "time",
    name_divergence = "divergence_from_baseline",
    name_timedifference = "time_from_baseline",
    assay.type="relabundance",
    FUN = vegan::vegdist,
    method="bray")


debug(.get_baseline_divergence)

tse2 <- addBaselineDivergence(
    tse,
    group = "bmi_group",
    time_field = "time",
    name_divergence = "divergence_from_baseline",
    name_timedifference = "time_from_baseline",
    assay.type="relabundance",
    FUN = vegan::vegdist,
    method="bray")

tse2 <- getBaselineDivergence(
    tse,
    baseline_sample = "Sample-875",
    group = "subject",
    time_field = "time",
    name_divergence = "divergence_from_baseline",
    name_timedifference = "time_from_baseline",
    assay.type="relabundance",
    FUN = vegan::vegdist,
    method="bray")

################################################################################

#library(miaTime)
library(TreeSummarizedExperiment)

data(hitchip1006)
tse <- mia::transformCounts(hitchip1006, method = "relabundance")

# Subset to speed up example
tse <- tse[, colData(tse)$subject %in% c("900", "934", "843", "875")]

# Using vegdist for divergence calculation, one can pass
# the dissimilarity method from the vegan::vegdist options
# via the "method" argument
debug(.get_stepwise_divergence2)
tse2 <- getStepwiseDivergence(tse, group = "subject",
                              time_interval = 1,
                              time_field = "time",
                              assay.type="relabundance",
                              FUN = vegan::vegdist,
                              method="bray")


####################
tse2 <- addStepwiseDivergence(tse, group = "subject",
                            time_interval = 1,
                                                            time_field = "time",
                                                            assay.type="relabundance",
                                                            FUN = vegan::vegdist,
                                                            method="bray")

