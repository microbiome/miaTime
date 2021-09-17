
data(GlobalPatterns, package="mia")
se <- GlobalPatterns

test_that("addRank", {
    expect_error(miaTime:::addRank(x = se, field = "SampleType", na.last = TRUE , ties.method = "first"),
                'argument "rank_field_name" is missing')

    se$newfield <- matrix(1:26, nrow = 26, ncol = 1)
    se2 <- addRank(se, field = "newfield" ,
                rank_field_name = "newfield_rank" ,
                na.last = TRUE , ties.method = "first")

    expect_equal(dim(se$newfield), c(26,1))
    expect_s4_class(se2, "TreeSummarizedExperiment")

})
