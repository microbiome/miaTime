# Check that input checks work
test_that("getSurvival validates input", {
    # Load and prepare data
    data(crohn_survival)
    tse <- crohn_survival
    tse <- transformAssay(tse, method = "relabundance")
    
    # Test that function works with valid inputs
    expect_no_error({
        fit <- getSurvival(
            tse, 
            assay.type = "relabundance",
            time.col = "event_time", 
            event.col = "event"
        )
    })
    
    expect_error({
        getSurvival(
            tse, 
            assay.type = "relabundance"
        )
    })
    expect_error({
        getSurvival(
            tse
        )
    })
    
    # Test error with non-existent time column
    expect_error({
        getSurvival(
            tse, 
            assay.type = "relabundance",
            time.col = "nonexistent_time", 
            event.col = "event"
        )
    })
    
    expect_error({
        getSurvival(
            tse, 
            assay.type = "relabundance",
            event.col = "event"
        )
    })
    
    # Test error with non-existent event column
    expect_error({
        getSurvival(
            tse, 
            assay.type = "relabundance",
            time.col = "event_time", 
            event.col = "nonexistent_event"
        )
    })
    
    expect_error({
        getSurvival(
            tse, 
            assay.type = "relabundance",
            time.col = "event_time"
        )
    })
    
    # Test error with non-existent assay
    expect_error({
        getSurvival(
            tse, 
            assay.type = "nonexistent_assay",
            time.col = "event_time", 
            event.col = "event"
        )
    })
})

# Check that method works correctly
test_that("getSurvival works correctly", {
    # Load the example data
    data(crohn_survival)
    tse <- crohn_survival
    
    # Check that the data loaded correctly
    expect_s4_class(tse, "TreeSummarizedExperiment")
    
    # Verify required columns exist
    expect_true("event_time" %in% colnames(colData(tse)))
    expect_true("event" %in% colnames(colData(tse)))
    
    # Transform to relative abundance as in example
    tse <- transformAssay(tse, method = "relabundance")
    
    # Test basic getSurvival functionality
    fit <- getSurvival(
        tse, 
        assay.type = "relabundance",
        time.col = "event_time", 
        event.col = "event"
    )
    
    # Check that result is a list with expected components
    expect_type(fit, "list")
    expect_true(all(c("coef", "risk_scores", "c_index", "fit") %in% names(fit)))
    
    # Check coefficients
    expect_type(fit$coef, "double")
    expect_true(is.numeric(fit$coef))
    
    # Check risk scores match number of samples
    expect_length(fit$risk_scores, ncol(tse))
    expect_true(is.numeric(fit$risk_scores))
    expect_true(all(is.finite(fit$risk_scores)))
    
    # Check C-index is valid
    expect_true(is.numeric(fit$c_index))
    expect_true(fit$c_index >= 0 && fit$c_index <= 1)
    
    # Check that fit object exists and is correct type
    expect_true(!is.null(fit$fit))
    expect_s3_class(fit$fit, "cv.glmnet")
    
    # Check cross-validation results are present (penalized = TRUE by default)
    expect_true("c_index_cv_mean" %in% names(fit))
    expect_true("c_index_cv_sd" %in% names(fit))
    expect_true(is.numeric(fit$c_index_cv_mean))
    expect_true(is.numeric(fit$c_index_cv_sd))
})

# Check that method handles different params passed in ...
test_that("getSurvival handles different parameters", {
    # Load and prepare data
    data(crohn_survival)
    tse <- crohn_survival
    tse <- transformAssay(tse, method = "relabundance")
    
    # Test with different alpha values
    fit_ridge <- getSurvival(
        tse, 
        assay.type = "relabundance",
        time.col = "event_time", 
        event.col = "event",
        alpha = 0.1  # Ridge regression
    )
    
    fit_lasso <- getSurvival(
        tse, 
        assay.type = "relabundance",
        time.col = "event_time", 
        event.col = "event", 
        alpha = 1.0  # Lasso regression
    )
    
    # Both should return valid results
    expect_type(fit_ridge, "list")
    expect_type(fit_lasso, "list")
    
    # Results should be different due to different regularization
    expect_false(identical(fit_ridge$coef, fit_lasso$coef))
    expect_false(identical(fit_ridge$risk_scores, fit_lasso$risk_scores))
    
    # Test with standard Cox regression (no penalization)
    fit_standard <- getSurvival(
        tse, 
        assay.type = "relabundance",
        time.col = "event_time", 
        event.col = "event",
        penalized = FALSE
    )
    
    expect_type(fit_standard, "list")
    expect_s3_class(fit_standard$fit, "coxph")
    
    # Standard model should not have CV results
    expect_false("c_index_cv_mean" %in% names(fit_standard))
    expect_false("c_index_cv_sd" %in% names(fit_standard))
})

# Check that method handles covariates
test_that("getSurvival handles covariates", {
    # Load and prepare data
    data(crohn_survival)
    tse <- crohn_survival
    tse <- transformAssay(tse, method = "relabundance")
    
    set.seed(42)  # For reproducibility
    
    # Determine sample size
    n <- ncol(tse)
    
    # Generate synthetic covariates
    colData(tse)$age <- round(rnorm(n, mean = 40, sd = 12))  # ages, e.g. 20–70
    colData(tse)$sex <- factor(sample(c("M", "F"), n, replace = TRUE))
    
    # Test with available covariates (use first available column)
    fit_with_covs <- getSurvival(
        tse, 
        assay.type = "relabundance",
        time.col = "event_time", 
        event.col = "event",
        col.var = c("age", "sex")
    )
    
    # Test without covariates for comparison
    fit_without_covs <- getSurvival(
        tse, 
        assay.type = "relabundance",
        time.col = "event_time", 
        event.col = "event"
    )
    
    # Both should return valid results
    expect_type(fit_with_covs, "list")
    expect_type(fit_without_covs, "list")
    
    # Results should be different when covariates are included
    expect_false(identical(fit_with_covs$risk_scores, fit_without_covs$risk_scores))
    
    # Both should have same structure
    essential_components <- c("coef", "risk_scores", "c_index", "fit")
    expect_true(all(essential_components %in% names(fit_with_covs)))
    expect_true(all(essential_components %in% names(fit_without_covs)))
})

# Check that method handles coefficient filtering
test_that("getSurvival handles coefficient filtering", {
    # Load and prepare data
    data(crohn_survival)
    tse <- crohn_survival
    tse <- transformAssay(tse, method = "relabundance")
    
    # Test with different coefficient thresholds
    fit_no_filter <- getSurvival(
        tse, 
        assay.type = "relabundance",
        time.col = "event_time", 
        event.col = "event",
        coef.threshold = 0  # No filtering
    )
    
    fit_filtered <- getSurvival(
        tse, 
        assay.type = "relabundance",
        time.col = "event_time", 
        event.col = "event",
        coef.threshold = 0.1  # Filter small coefficients
    )
    
    # Filtered result should have fewer or equal coefficients
    expect_true(length(fit_filtered$coef) <= length(fit_no_filter$coef))
    
    # All remaining coefficients should exceed the threshold
    if(length(fit_filtered$coef) > 0) {
        expect_true(all(abs(fit_filtered$coef) > 0.1))
    }
    
    # Test with high threshold that might filter out all coefficients
    fit_high_filter <- getSurvival(
        tse, 
        assay.type = "relabundance",
        time.col = "event_time", 
        event.col = "event",
        coef.threshold = 10  # Very high threshold
    )
    
    # Should still return valid structure even if no coefficients pass threshold
    expect_type(fit_high_filter, "list")
    expect_true("coef" %in% names(fit_high_filter))
})

# Check that method handles different transformations
test_that("getSurvival handles different assay types", {
    # Load data
    data(crohn_survival)
    tse <- crohn_survival
    
    # Test with original counts
    fit_counts <- getSurvival(
        tse, 
        assay.type = "counts",
        time.col = "event_time", 
        event.col = "event"
    )
    
    # Transform and test with relative abundances
    tse <- transformAssay(tse, method = "relabundance")
    fit_relabundance <- getSurvival(
        tse, 
        assay.type = "relabundance",
        time.col = "event_time", 
        event.col = "event"
    )
    
    # Both should return valid results
    expect_type(fit_counts, "list")
    expect_type(fit_relabundance, "list")
    
    # Results should be different due to different data scaling
    expect_false(identical(fit_counts$coef, fit_relabundance$coef))
    expect_false(identical(fit_counts$risk_scores, fit_relabundance$risk_scores))
    
    # Add CLR transformation and test
    tse <- transformAssay(tse, method = "clr", pseudocount = TRUE)
    fit_clr <- getSurvival(
        tse, 
        assay.type = "clr",
        time.col = "event_time", 
        event.col = "event"
    )
    
    expect_type(fit_clr, "list")
    expect_false(identical(fit_clr$coef, fit_relabundance$coef))
})

# Check that methods results are reproducible
test_that("getSurvival results are reproducible", {
    # Load and prepare data
    data(crohn_survival)
    tse <- crohn_survival
    tse <- transformAssay(tse, method = "relabundance")
    
    # Set seed and run analysis
    set.seed(123)
    fit1 <- getSurvival(
        tse, 
        assay.type = "relabundance",
        time.col = "event_time", 
        event.col = "event",
        nfolds = 5  # Use fewer folds for faster testing
    )
    
    # Set same seed and run analysis again
    set.seed(123)
    fit2 <- getSurvival(
        tse, 
        assay.type = "relabundance",
        time.col = "event_time", 
        event.col = "event",
        nfolds = 5
    )
    
    # Results should be identical with same seed
    expect_identical(fit1$coef, fit2$coef)
    expect_identical(fit1$risk_scores, fit2$risk_scores)
    expect_equal(fit1$c_index, fit2$c_index)
})