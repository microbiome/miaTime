test_that("getVariogramModel input checks", {
    data(hitchip1006, package = "miaTime")
    tse <- hitchip1006
    expect_error(getVariogramModel(tse, time.col = "x", group = "subject"))
    expect_error(getVariogramModel(tse, time.col = "time", group = "x"))
    expect_error(getVariogramModel(tse, time.col = "time",
        group = "subject", form = "spherical"))
    expect_error(getVariogramModel(tse, time.col = "time",
        group = "subject", conf.level = 0))
})

# Pairs of 40 subjects whose dissimilarities follow an exponential variogram
# with nugget 0.1, partial sill 0.5 and timescale 10, with gamma noise
.make_exponential_pairs <- function(){
    set.seed(1)
    h <- rep(c(0, 1, 2, 4, 7, 10, 15, 20, 30, 40), 40)
    g <- 0.1 + 0.5 * (1 - exp(-h / 10))
    data.frame(subject = rep(paste0("s", 1:40), each = 10), h = h, mid = 20,
        d = g * rgamma(400, shape = 50, rate = 50))
}

test_that("The exponential form is recovered", {
    pairs <- .make_exponential_pairs()
    fit <- miaTime:::.fit_variogram(pairs, "exponential")
    expect_true(fit$converged)
    res <- miaTime:::.variogram_parameters(fit, 0.95)
    est <- setNames(res$estimate, res$parameter)
    expect_equal(unname(est[c("nugget", "partial.sill", "sill", "timescale")]),
        c(0.1, 0.5, 0.6, 10), tolerance = 0.1)
    expect_true(all(res$lower < res$estimate & res$estimate < res$upper))
    expect_true(res$identified[res$parameter == "timescale"])
    expect_true(all(is.na(res$identified[res$parameter != "timescale"])))
    expect_equal(res$df[[1]], 39)
    # The simplest adequate form is chosen
    set.seed(1)
    cv <- miaTime:::.cv_forms(pairs, c("exponential", "power", "linear"))
    expect_identical(cv$form[cv$chosen], "exponential")
})

test_that("getVariogramModel output", {
    data(hitchip1006, package = "miaTime")
    tse <- transformAssay(hitchip1006, method = "hellinger")
    args <- list(tse, time.col = "time", group = "subject",
        assay.type = "hellinger", method = "euclidean", squared = TRUE)
    res <- do.call(getVariogramModel, c(args, form = "exponential"))
    expect_s4_class(res, "DataFrame")
    expect_identical(res$parameter,
        c("nugget", "partial.sill", "sill", "timescale"))
    expect_true(all(res$form == "exponential"))
    expect_null(metadata(res)$cv)
    expect_equal(res$estimate[3], res$estimate[1] + res$estimate[2])
    expect_identical(do.call(getVariogramModel, c(args, form = "linear"))$parameter,
        c("nugget", "slope"))
    # Cross-validation
    set.seed(1)
    res <- do.call(getVariogramModel, c(args,
        list(form = c("exponential", "linear"))))
    cv <- metadata(res)$cv
    expect_identical(cv$form, c("exponential", "linear"))
    expect_identical(sum(cv$chosen), 1L)
    expect_identical(res$form[[1]], cv$form[cv$chosen])
    # Stored in metadata
    tse <- do.call(addVariogramModel, c(args, form = "exponential",
        name = "model"))
    expect_equal(metadata(tse)[["model"]]$estimate,
        do.call(getVariogramModel, c(args, form = "exponential"))$estimate)
})
