library(testthat)
library(miaTime)

data(airway, package="airway")
se <- airway
newmatrix <- matrix(1:8, nrow = 8, ncol = 1)

test_that("addRank", {
    expect_error(miaTime:::addRank(x = se, field = "SampleName"),
                'argument "rank_field_name" is missing' , fixed = TRUE)
    expect_error(miaTime:::addRank(x = se, rank_field_name = "SampleName_rank"),
                'argument "field" is missing')

    Added_field <- addRank(se, field = "NewMatrixAdded" ,
                field_matrix = newmatrix, rank_field_name = "newRankedField" ,
                na.last = TRUE , ties.method = "first")

    expect_true(is.matrix(newmatrix))
    expect_equal(dim(newmatrix), c(8,1))
    expect_s4_class(Added_field, "RangedSummarizedExperiment")

})
