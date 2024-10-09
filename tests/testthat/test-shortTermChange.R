test_that("shortTermChange", {
    library(SummarizedExperiment)
    # Load dataset
    data(minimalgut)
    tse <- minimalgut
    # Check if the function handles empty input
    empty_se <- SummarizedExperiment()
    expect_error(shortTermChange(empty_se), 
                 "No data available in `x`")
    # Check if assay.type argument works
    # tse_invalid <- tse
    # expect_error(
    #     shortTermChange(tse_invalid, assay.type = "invalid_assay"),
    #     "'assay.type' must be a valid name of assays(x)"
    # )
    # Check that rarefy and compositional cannot both be TRUE
    expect_error(shortTermChange(tse, rarefy = TRUE, compositional = TRUE), 
                 "Both rarefy and compositional cannot be TRUE simultaneously")
    # Check if the depth argument is greater than minimum counts
    min_depth <- min(assay(tse, "counts"))
    expect_error(shortTermChange(tse, depth = min_depth + 1),
                 "Depth cannot be greater than the minimum number of counts in your data")
    # Check if rarefy = TRUE works
    result <- shortTermChange(tse, rarefy = TRUE)
    expect_true(is.data.frame(result))
    # Check if compositional = TRUE works
    result <- shortTermChange(tse, compositional = TRUE)
    expect_true(is.data.frame(result))  
    # Expected output should be a ggplot object
    result <- shortTermChange(tse, plot = TRUE)
    expect_s3_class(result, "gg")  
    # Should still return a dataframe
    result <- shortTermChange(tse, rarefy = TRUE, 
                              compositional = FALSE, 
                              additional_arg = "value")
    expect_true(is.data.frame(result))  
    
    short_time_labels <- c("74.5h", "173h", "438h", "434h", "390h")
    # Subset samples by Time_label and StudyIdentifier 
    tse_filtered <- tse[, !(colData(tse)$Time_label %in% short_time_labels)]
    tse_filtered <- tse_filtered[, (colData(tse_filtered)$StudyIdentifier == "Bioreactor A")]
    
    expect_true(all(!(colData(tse_filtered)$Time_label %in% short_time_labels)))
    # Expected output should be a ggplot object
    result <- shortTermChange(tse_filtered, 
                              rarefy = TRUE, plot = TRUE)
    expect_s3_class(result, "gg")  
    
    result <- shortTermChange(tse_filtered, 
                              rarefy = FALSE, plot = FALSE)
    # Expected output is a dataframe
    expect_true(is.data.frame(result))  
    # Ensure dataframe has expected columns
    expect_true("OTU" %in% colnames(result))  
    expect_true("growth_diff" %in% colnames(result))
    
    result <- shortTermChange(tse_filtered, 
                              rarefy = TRUE, plot = FALSE)
    
    result <- shortTermChange(tse_filtered, 
                              compositional = TRUE, plot = FALSE)
    # Expected output is a dataframe
    expect_true(is.data.frame(result))  
    
    result <- shortTermChange(tse_filtered, rarefy = FALSE, 
                              compositional = FALSE, 
                              plot = FALSE)
    expect_true("growth_diff" %in% colnames(result))
    # Test some expected properties (e.g., that growth_diff isn't all NAs)
    expect_false(all(is.na(result$growth_diff)))
    
    min_depth <- min(assay(tse_filtered, "counts"))
    result <- shortTermChange(tse_filtered, rarefy = TRUE, depth = min_depth)
    expect_true(is.data.frame(result)) 
    expect_error(shortTermChange(tse_filtered, 
                                 rarefy = TRUE, 
                                 depth = min_depth + 1),
                 "Depth cannot be greater than the minimum number of counts in your data")
    
})