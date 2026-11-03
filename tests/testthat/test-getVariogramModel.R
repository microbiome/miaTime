test_that("getVariogramModel input checks", {
    data(hitchip1006, package = "miaTime")
    tse <- hitchip1006
    expect_error(getVariogramModel(tse, time.col = "x", group = "subject"))
    expect_error(getVariogramModel(tse, time.col = "time", group = "x"))
    expect_error(getVariogramModel(tse, time.col = "time",
        group = "subject", form = "spherical"))
    expect_error(getVariogramModel(tse, time.col = "time",
        group = "subject", conf.level = 0))
    expect_error(getVariogramModel(tse, time.col = "time",
        group = "subject", split.by = "sex", bootstrap = 0))
    expect_error(getVariogramModel(tse, time.col = "time",
        group = "subject", split.by = "sex", BPPARAM = 1))
    expect_error(getVariogramModel(tse, time.col = "time",
        group = "subject", time.point = "median"))
    tse$form <- tse$sex
    expect_error(getVariogramModel(tse, time.col = "time",
        group = "subject", split.by = "form"), "output column")
    # The model compares subject-level groups only
    tse$state <- rep(c("a", "b"), length.out = ncol(tse))
    expect_error(suppressMessages(getVariogramModel(tse, time.col = "time",
        group = "subject", split.by = "state")))
})

# Pairs of 40 subjects whose dissimilarities follow an exponential variogram
# with partial sill 0.5 and timescale 10, with gamma noise
.make_exponential_pairs <- function(nugget = 0.1, prefix = "s"){
    h <- rep(c(0, 1, 2, 4, 7, 10, 15, 20, 30, 40), 40)
    g <- nugget + 0.5 * (1 - exp(-h / 10))
    data.frame(subject = rep(paste0(prefix, 1:40), each = 10), h = h,
        mid = 20, d = g * rgamma(400, shape = 50, rate = 50))
}

test_that("The exponential form is recovered", {
    set.seed(1)
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

test_that("Parameters are compared between groups", {
    # The second group has three times the nugget
    set.seed(1)
    pairs <- rbind(
        cbind(.make_exponential_pairs(0.1, "a"), level = "A"),
        cbind(.make_exponential_pairs(0.3, "b"), level = "B"))
    pairs$level <- factor(pairs$level)
    res <- miaTime:::.compare_variogram_models(pairs, "exponential", 0.95,
        99, BiocParallel::SerialParam())
    ratio <- res[res$type == "ratio", ]
    expect_identical(unique(res$level), c("A", "B", "B / A"))
    nug <- ratio[ratio$parameter == "nugget", ]
    expect_equal(nug$estimate, 3, tolerance = 0.2)
    expect_true(nug$lower > 1 && nug$p.value < 0.05)
    ps <- ratio[ratio$parameter == "partial.sill", ]
    expect_true(ps$lower < 1 && ps$upper > 1)
    expect_true(all(res$lower[res$type == "parameter"] <=
        res$estimate[res$type == "parameter"]))
    expect_true(all(res$max.lag == 40))
    # Timescales are compared only when identified in both groups
    ts <- res[res$parameter == "timescale", ]
    expect_identical(is.na(ts$p.value[3]),
        !(ts$identified[1] && ts$identified[2]))

    data(hitchip1006, package = "miaTime")
    tse <- transformAssay(hitchip1006, method = "hellinger")
    set.seed(1)
    res <- getVariogramModel(tse, time.col = "time", group = "subject",
        assay.type = "hellinger", method = "euclidean", squared = TRUE,
        form = "linear", split.by = "sex", bootstrap = 19, time.point = 5)
    expect_identical(unique(res$sex), c("female", "male", "male / female"))
    expect_identical(res$type, rep(c("parameter", "ratio"), c(4, 2)))
    expect_true(all(res$time.point == 5))
    expect_true(all(is.na(res$p.value[res$type == "parameter"])))
})

test_that("Permutation test for temporal structure", {
    # Twenty subjects with eight samples each: a random walk has temporal
    # structure, independent values have none
    time <- rep(0:7, 20)
    subject <- rep(paste0("s", 1:20), each = 8)
    ij <- do.call(rbind, lapply(split(seq_along(time), subject), function(i){
        t(combn(i, 2))
    }))
    .pairs <- function(x){
        data.frame(subject = subject[ij[, 1]],
            h = abs(time[ij[, 2]] - time[ij[, 1]]),
            mid = (time[ij[, 1]] + time[ij[, 2]]) / 2,
            d = (x[ij[, 1]] - x[ij[, 2]])^2, k = ij[, 1], l = ij[, 2])
    }
    set.seed(1)
    walk <- ave(rnorm(160), subject, FUN = cumsum)
    res <- miaTime:::.temporal_structure(.pairs(walk), time, 99)
    expect_true(res$statistic > 0)
    expect_equal(res$p.value, 0.01)
    res <- miaTime:::.temporal_structure(.pairs(rnorm(160)), time, 99)
    expect_true(res$p.value > 0.05)

    data(hitchip1006, package = "miaTime")
    tse <- transformAssay(hitchip1006, method = "hellinger")
    args <- list(tse, time.col = "time", group = "subject",
        assay.type = "hellinger", method = "euclidean", squared = TRUE,
        form = "linear")
    expect_false("test" %in% do.call(getVariogramModel, args)$type)
    set.seed(1)
    res <- do.call(getVariogramModel, c(args, permutations = 19))
    tt <- res[res$type == "test", ]
    expect_identical(tt$parameter, "temporal.structure")
    expect_true(tt$p.value >= 1 / 20 && tt$p.value <= 1)
    set.seed(1)
    res <- do.call(getVariogramModel, c(args, permutations = 19,
        split.by = "sex", bootstrap = 9))
    expect_identical(res$sex[res$type == "test"], c("female", "male"))
    expect_error(do.call(getVariogramModel, c(args, permutations = -1)))
})
