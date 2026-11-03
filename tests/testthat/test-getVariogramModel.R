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
    fit <- miaTime:::.fit_variogram(pairs, miaTime:::.form_spec("exponential"))
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
    cv <- miaTime:::.cv_forms(pairs, lapply(c(exponential = "exponential",
        power = "power", linear = "linear"), miaTime:::.form_spec))
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
        c("nugget", "partial.sill", "sill", "timescale", "stationarity"))
    expect_true(all(res$form == "exponential"))
    expect_null(metadata(res)$cv)
    expect_equal(res$estimate[3], res$estimate[1] + res$estimate[2])
    expect_identical(do.call(getVariogramModel, c(args, form = "linear"))$parameter,
        c("nugget", "slope", "stationarity"))
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
    res <- miaTime:::.compare_variogram_models(pairs,
        list(exponential = miaTime:::.form_spec("exponential")), 0.95,
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
    expect_identical(res$type, rep(c("parameter", "ratio", "test"),
        c(4, 2, 2)))
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
    expect_false("temporal.structure" %in%
        do.call(getVariogramModel, args)$parameter)
    set.seed(1)
    res <- do.call(getVariogramModel, c(args, permutations = 19))
    tt <- res[res$parameter == "temporal.structure", ]
    expect_identical(tt$parameter, "temporal.structure")
    expect_true(tt$p.value >= 1 / 20 && tt$p.value <= 1)
    set.seed(1)
    res <- do.call(getVariogramModel, c(args, permutations = 19,
        split.by = "sex", bootstrap = 9))
    expect_identical(res$sex[res$parameter == "temporal.structure"],
        c("female", "male"))
    expect_error(do.call(getVariogramModel, c(args, permutations = -1)))
})

test_that("Stationarity check", {
    # Relative residuals that grow by 0.01 per unit of midpoint time, with
    # pairs of every subject spread over time
    set.seed(1)
    pairs <- .make_exponential_pairs()
    pairs$mid <- rep(seq(0, 90, by = 10), 40)
    fit <- miaTime:::.fit_variogram(pairs, miaTime:::.form_spec("exponential"))
    res <- miaTime:::.stationarity(fit, 0.95)
    expect_true(abs(res$estimate) < 3 * res$se)
    fit$pairs$d <- fit$fitted * (1 + 0.01 * (pairs$mid - 45) +
        rnorm(400, sd = 0.05))
    res <- miaTime:::.stationarity(fit, 0.95)
    expect_equal(res$estimate, 0.01, tolerance = 0.05)
    expect_true(res$p.value < 0.001)
    expect_true(res$lower < 0.01 && res$upper > 0.01)
    expect_identical(res$df, 39)
    # No slope when every pair has the same midpoint
    fit$pairs$mid <- 20
    expect_true(is.na(miaTime:::.stationarity(fit, 0.95)$estimate))
})

test_that("Seasonal forms add an annual cycle", {
    # Pairs whose dissimilarity rises and falls with a period of 365 days
    set.seed(1)
    h <- rep(seq(0, 1460, by = 20), 30)
    g <- 0.1 + 0.0002 * h + 0.2 * (1 - cos(2 * pi * h / 365))
    pairs <- data.frame(subject = rep(paste0("s", 1:30), each = 74), h = h,
        mid = 800, d = g * rgamma(length(h), shape = 50, rate = 50))
    fit <- miaTime:::.fit_variogram(pairs,
        miaTime:::.form_spec("linear.seasonal", period = 365))
    res <- miaTime:::.variogram_parameters(fit, 0.95)
    est <- setNames(res$estimate, res$parameter)
    expect_equal(unname(est[c("nugget", "slope", "seasonal")]),
        c(0.1, 0.0002, 0.2), tolerance = 0.1)
    set.seed(1)
    forms <- lapply(c(linear = "linear",
        linear.seasonal = "linear.seasonal"), miaTime:::.form_spec,
        period = 365)
    cv <- miaTime:::.cv_forms(pairs, forms)
    expect_identical(cv$form[cv$chosen], "linear.seasonal")
    data(hitchip1006, package = "miaTime")
    expect_error(getVariogramModel(hitchip1006, time.col = "time",
        group = "subject", form = "linear.seasonal"), "period")
    expect_error(getVariogramModel(hitchip1006, time.col = "time",
        group = "subject", form = "linear.seasonal", period = -1))
})

test_that("Technical replicates split the nugget", {
    # Forty subjects; technical noise 0.05, fast change 0.1, partial sill 0.5
    # and timescale 10. Replicate pairs (lag 0) carry the technical part only
    set.seed(1)
    pairs <- .make_exponential_pairs(nugget = 0.15)
    pairs <- pairs[pairs$h > 0, ]
    pairs$rep <- FALSE
    reps <- data.frame(subject = paste0("s", 1:40), h = 0, mid = 20,
        d = 0.05 * rgamma(40, shape = 50, rate = 50), rep = TRUE)
    pairs <- rbind(pairs, reps)
    fit <- miaTime:::.fit_variogram(pairs,
        miaTime:::.form_spec("exponential", replicates = TRUE))
    res <- miaTime:::.variogram_parameters(fit, 0.95)
    est <- setNames(res$estimate, res$parameter)
    expect_equal(unname(est[c("nugget.technical", "nugget.fast", "nugget")]),
        c(0.05, 0.1, 0.15), tolerance = 0.1)
    # Replicate pairs are not used by the stationarity check or the lag range
    expect_identical(miaTime:::.stationarity(fit, 0.95)$n.pairs,
        sum(!pairs$rep))
    expect_identical(res$min.lag[[1]], 1)

    # From a SummarizedExperiment: two replicates of the first sample of
    # every subject
    set.seed(1)
    n <- 160
    mat <- matrix(runif(5 * n), 5, dimnames = list(paste0("f", 1:5), NULL))
    cd <- DataFrame(subject = rep(1:20, each = 8),
        time = rep(0:7, 20) + c(0, 0, runif(6, 0.1, 0.5)),
        rg = rep(c("r", "r", rep(NA, 6)), 20))
    cd$rg[!is.na(cd$rg)] <- paste0(cd$rg[!is.na(cd$rg)], cd$subject[
        !is.na(cd$rg)])
    tse <- SummarizedExperiment(assays = list(x = mat), colData = cd)
    res <- getVariogramModel(tse, time.col = "time", group = "subject",
        assay.type = "x", method = "euclidean", squared = TRUE,
        form = "linear", replicate.by = "rg")
    expect_true(all(c("nugget.technical", "nugget.fast") %in% res$parameter))
    expect_identical(res$n.pairs[[1]], 20L * 21L + 20L)
    expect_error(getVariogramModel(tse, time.col = "time",
        group = "subject", assay.type = "x", replicate.by = "x"))
})
