test_that("getTsMeasure", {
    data(minimalgut)
    tse <- minimalgut
    
    # Calculate time series measures
    result <- getTsMeasure(tse, time.field = "Time.hr", measures = c("acf", "pacf", 
       "Box.test","arima"))
    
    # Test that measures are calculated and correctly stored in the list
    expect_true("acf" %in% names(result))
    expect_true("pacf" %in% names(result))
    expect_true("Box.test" %in% names(result))
    expect_true("arima" %in% names(result))
    expect_type(result$acf, "double")
    expect_type(result$pacf, "double")
    expect_s3_class(result$Box.test, "htest")
    expect_s3_class(result$arima, "Arima")
}
)