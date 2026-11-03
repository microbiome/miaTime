#' @name
#' getVariogramModel
#'
#' @export
#'
#' @title
#' Nugget, sill and timescale of community change (variogram model)
#'
#' @description
#' Fits a form to the temporal variogram \eqn{g(h)}, the expected
#' dissimilarity between two samples of the same subject taken \eqn{h} time
#' units apart, and reports its nugget, sill and timescale. The nugget
#' separates measurement noise and change faster than the shortest lag from
#' slower change.
#'
#' @details
#' The forms are
#' \describe{
#'   \item{exponential}{\eqn{g(h) = c_0 + c_1 (1 - e^{-h / \phi})}: change
#'   that levels off at the sill \eqn{c_0 + c_1}, with timescale \eqn{\phi}.}
#'   \item{two.timescale}{\eqn{g(h) = c_0 + c_1 (1 - e^{-h / \phi_1}) +
#'   c_2 (1 - e^{-h / \phi_2})}, with \eqn{\phi_2 > \phi_1}.}
#'   \item{power}{\eqn{g(h) = c_0 + b h^\alpha}, \eqn{0 < \alpha \le 2}:
#'   change without a sill.}
#'   \item{linear}{\eqn{g(h) = c_0 + b h}.}
#' }
#' In every form \eqn{c_0 = g(0)} is the nugget: the dissimilarity between
#' samples taken at the same time, from measurement noise and change faster
#' than the shortest lag. Pairs of samples taken at the same time, such as
#' technical replicates, inform it directly.
#'
#' The form is fitted to all pairs of samples of the same subject by weighted
#' least squares with variance proportional to \eqn{g(h)^2}, with iterated
#' weights (Cressie 1985); the pairs of each subject are weighted so that
#' every subject counts equally. Pairs with lags longer than the longest
#' lag reached by at least \eqn{\max(5, 20\%)} of the subjects are left out,
#' because they come from few subjects and are dominated by them.
#'
#' When several forms are given, the form is chosen by ten-fold
#' cross-validation over subjects: each form is fitted to the subjects of
#' nine folds, and the held-out subjects' variogram, binned at quantiles of
#' the lags, is compared with the fitted one by the mean squared relative
#' difference. The simplest form whose loss is within one standard error of
#' the smallest is chosen. The losses are in \code{metadata(res)$cv}. Set a
#' seed (\code{set.seed()}) for reproducible folds.
#'
#' Standard errors are cluster-robust by subject (CR1) and are computed on
#' the log scale of each parameter, by the delta method; intervals use a t
#' distribution with \eqn{G - 1} degrees of freedom, where \eqn{G} is the
#' number of subjects. A timescale is identified only when its interval lies
#' within the fitted lag range, from the shortest positive lag to the
#' longest lag used; otherwise the data do not determine it, and change at
#' stated lags
#' (\code{\link[miaTime:getVariogram]{getVariogram()}}) is the quantity to
#' report.
#'
#' \strong{Groups.} Two groups can differ in change at a lag because of
#' measurement noise and fast change (the nugget) or because of slower
#' dynamics (the sill and the timescale). With \code{split.by}, the form is
#' fitted separately to the subjects of each level of a subject-level
#' variable, over the same lag range (the shortest of the levels' ranges),
#' and each level is compared with the first by the ratio of each parameter.
#' When several forms are given, the form is chosen once, by cross-validation
#' on the pairs of all levels, so that the levels share it. Intervals are
#' percentile intervals from \code{bootstrap} resamples of the subjects
#' within each level, each refitted; resamples whose fit fails or does not
#' converge are left out. The p-value of a ratio is twice the smaller share
#' of resamples on either side of 1. Timescales are compared only when they
#' are identified in both levels. P-values are not adjusted for multiple
#' comparisons. Set a seed (\code{set.seed()}) for reproducible intervals.
#'
#' \strong{Developing communities.} When change depends on the time since a
#' common origin (age in infant cohorts), groups sampled at different ages
#' differ even with the same development. With \code{time.point}, only pairs
#' whose midpoint \eqn{m} satisfies \eqn{(1 + m) / (1 + a)} in
#' \eqn{[1/2, 2]} are fitted, as in
#' \code{\link[miaTime:getVariogram]{getVariogram()}}, so that the levels
#' are compared around the same time \eqn{a}.
#'
#' The dissimilarities are those of
#' \code{\link[miaTime:getVariogram]{getVariogram()}}; see its details for
#' the choice of dissimilarity and for bounded dissimilarities.
#'
#' @references
#' Cressie N. (1985) Fitting variogram models by weighted least squares.
#' Mathematical Geology 17, 563–586. doi: 10.1007/BF01032109
#'
#' Diggle P.J. (1988) An approach to the analysis of repeated measurements.
#' Biometrics 44, 959–971. doi: 10.2307/2531727
#'
#' @return
#' \code{getVariogramModel} returns a \code{DataFrame} with one row per
#' parameter (and level of \code{split.by}, and comparison) and the columns
#' \describe{
#'   \item{\code{split.by}}{When \code{split.by} is given: the level, or the
#'   comparison (\code{"B / A"}).}
#'   \item{type}{\code{"parameter"}, or \code{"ratio"} for comparisons
#'   between levels.}
#'   \item{form}{The fitted form.}
#'   \item{parameter}{\code{"nugget"}, \code{"partial.sill"},
#'   \code{"sill"} and \code{"timescale"} (with suffixes \code{.1} and
#'   \code{.2} for the two-timescale form), \code{"scale"} and
#'   \code{"exponent"} (power form) or \code{"slope"} (linear form).}
#'   \item{estimate, lower, upper}{The estimate, or the ratio between
#'   levels, with its confidence interval.}
#'   \item{se}{Standard error on the log scale (without \code{split.by}).}
#'   \item{df}{Degrees of freedom of the t distribution (without
#'   \code{split.by}).}
#'   \item{p.value}{For comparisons, the p-value of the test that the ratio
#'   is 1.}
#'   \item{identified}{For timescales, whether the interval lies within the
#'   fitted lag range.}
#'   \item{n.subjects, n.pairs}{Subjects and pairs in the fit.}
#'   \item{min.lag, max.lag}{The fitted lag range.}
#'   \item{time.point}{The time around which the pairs were fitted, when
#'   \code{time.point} is given.}
#' }
#' \code{addVariogramModel} returns \code{x} with this table in its
#' \code{metadata}.
#'
#' @inheritParams getVariogram
#'
#' @param form \code{Character vector}. The forms among which one is chosen
#' by cross-validation: \code{"exponential"}, \code{"two.timescale"},
#' \code{"power"} and \code{"linear"}. See details. (Default: all four)
#'
#' @param split.by \code{NULL} or \code{Character scalar}. Specifies a name
#' of a column from \code{colData} with one value per subject. The form is
#' fitted to each of its levels, and each level is compared with the first.
#' (Default: \code{NULL})
#'
#' @param time.point \code{NULL} or \code{Numeric scalar}. The time around
#' which pairs are fitted, for communities whose change depends on time since
#' a common origin. See details. (Default: \code{NULL})
#'
#' @param bootstrap \code{Integer scalar}. The number of resamples of
#' subjects within each level of \code{split.by}. (Default: \code{999})
#'
#' @param BPPARAM A
#' \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#' object specifying whether the resamples are fitted in parallel.
#' (Default: \code{BiocParallel::SerialParam()})
#'
#' @param name \code{Character scalar}. Specifies a name for storing the
#' results in \code{metadata}. (Default: \code{"variogram_model"})
#'
#' @param ... additional arguments passed to
#' \code{\link[mia:getDissimilarity]{getDissimilarity()}}.
#'
#' @examples
#' library(miaTime)
#'
#' data(hitchip1006)
#' tse <- transformAssay(hitchip1006, method = "hellinger")
#'
#' # Nugget, sill and timescale of the exponential form, as squared
#' # Hellinger distance
#' res <- getVariogramModel(
#'     tse,
#'     time.col = "time",
#'     group = "subject",
#'     assay.type = "hellinger",
#'     method = "euclidean",
#'     squared = TRUE,
#'     form = "exponential"
#'     )
#' res
#'
#' # Choose among three forms by cross-validation and store the results in
#' # metadata
#' set.seed(1)
#' tse <- addVariogramModel(
#'     tse,
#'     time.col = "time",
#'     group = "subject",
#'     assay.type = "hellinger",
#'     method = "euclidean",
#'     squared = TRUE,
#'     form = c("exponential", "power", "linear")
#'     )
#' metadata(tse)[["variogram_model"]]
#'
#' # Compare the sexes; use more resamples in practice
#' set.seed(1)
#' res <- getVariogramModel(
#'     tse,
#'     time.col = "time",
#'     group = "subject",
#'     assay.type = "hellinger",
#'     method = "euclidean",
#'     squared = TRUE,
#'     form = "linear",
#'     split.by = "sex",
#'     bootstrap = 49
#'     )
#' res[res$type == "ratio", ]
#'
#' @seealso
#' \code{\link[miaTime:getVariogram]{getVariogram()}}
#'
NULL

#' @rdname getVariogramModel
#' @export
setMethod("addVariogramModel", signature = c(x = "SummarizedExperiment"),
    function(x, name = "variogram_model", ...){
        temp <- .check_input(name, list("character scalar"))
        x <- .check_and_get_altExp(x, ...)
        args <- c(list(x = x), list(...)[!names(list(...)) %in% c("altexp")])
        res <- do.call(getVariogramModel, args)
        x <- .add_values_to_metadata(x, name, res, ...)
        return(x)
    }
)

#' @rdname getVariogramModel
#' @export
setMethod("getVariogramModel", signature = c(x = "SummarizedExperiment"),
    function(
        x,
        time.col,
        group,
        assay.type = "counts",
        method = "bray",
        form = c("exponential", "two.timescale", "power", "linear"),
        split.by = NULL,
        time.point = NULL,
        squared = FALSE,
        bounded = FALSE,
        conf.level = 0.95,
        bootstrap = 999,
        BPPARAM = BiocParallel::SerialParam(),
        ...){
        ############################# INPUT CHECK ##############################
        x <- .check_and_get_altExp(x, ...)
        temp <- .check_input(
            time.col, list("character scalar"), colnames(colData(x)))
        if( !is.numeric(x[[time.col]]) ){
            stop("'time.col' must specify numeric column from colData(x)",
                call. = FALSE)
        }
        temp <- .check_input(
            group, list("character scalar"), colnames(colData(x)))
        .check_assay_present(assay.type, x)
        temp <- .check_input(method, list("character scalar"))
        temp <- .check_input(form, list("character vector"),
            names(.VARIOGRAM_FORMS))
        temp <- .check_input(
            split.by, list(NULL, "character scalar"), colnames(colData(x)))
        .check_split_by_name(split.by, c(.VARIOGRAM_COLUMNS, "form",
            "parameter", "identified", "min.lag", "max.lag"))
        temp <- .check_input(time.point, list(NULL, "numeric scalar"))
        temp <- .check_input(squared, list("logical scalar"))
        temp <- .check_input(bounded, list("logical scalar"))
        temp <- .check_input(conf.level, list("numeric scalar"))
        if( conf.level <= 0 || conf.level >= 1 ){
            stop("'conf.level' must be between 0 and 1.", call. = FALSE)
        }
        temp <- .check_input(bootstrap, list("numeric scalar"))
        if( bootstrap < 1 || bootstrap != round(bootstrap) ){
            stop("'bootstrap' must be a positive integer.", call. = FALSE)
        }
        if( !is(BPPARAM, "BiocParallelParam") ){
            stop("'BPPARAM' must be a BiocParallelParam object.",
                call. = FALSE)
        }
        ########################### INPUT CHECK END ############################
        args <- list(...)[!names(list(...)) %in% c("altexp")]
        pairs <- .get_within_subject_pairs(
            x, time.col, group, assay.type, method, split.by, squared,
            bounded, args)
        if( isTRUE(attr(pairs, "pair.level")) ){
            stop("'split.by' must have one value per subject.", call. = FALSE)
        }
        if( !is.null(time.point) ){
            r <- (1 + pmax(pairs$mid, 0)) / (1 + max(time.point, 0))
            pairs <- pairs[r >= 1 / 2 & r <= 2, , drop = FALSE]
        }
        form <- unique(form)
        if( is.null(split.by) ){
            pairs <- pairs[pairs$h <= .max_lag(pairs), , drop = FALSE]
            cv <- NULL
            if( length(form) > 1L ){
                cv <- .cv_forms(pairs, form)
                form <- cv$form[cv$chosen]
            }
            fit <- .fit_variogram(pairs, form)
            if( !fit$converged ){
                warning("The fit of the ", form, " form did not converge.",
                    call. = FALSE)
            }
            res <- .variogram_parameters(fit, conf.level)
        } else{
            res <- .compare_variogram_models(
                pairs, form, conf.level, bootstrap, BPPARAM)
            cv <- attr(res, "cv")
            names(res)[names(res) == "level"] <- split.by
        }
        if( !is.null(time.point) ){
            res$time.point <- time.point
        }
        res <- DataFrame(res, check.names = FALSE)
        metadata(res)$cv <- cv
        return(res)
    }
)

################################ HELP FUNCTIONS ################################

# Each form: names of its coefficients (positive parameters on the log scale,
# the exponent of the power form on the logit scale of alpha / 2), the
# variogram given the natural parameters, and the reported parameters
.VARIOGRAM_FORMS <- list(
    exponential = list(
        coef = c("log_c0", "log_c1", "log_phi"),
        g = function(h, P) P$c0 + P$c1 * (1 - exp(-h / P$phi)),
        report = function(P) c(
            nugget = P$c0, partial.sill = P$c1, sill = P$c0 + P$c1,
            timescale = P$phi)),
    two.timescale = list(
        coef = c("log_c0", "log_c1", "log_phi1", "log_c2", "log_r"),
        g = function(h, P) P$c0 + P$c1 * (1 - exp(-h / P$phi1)) +
            P$c2 * (1 - exp(-h / P$phi2)),
        report = function(P) c(
            nugget = P$c0, partial.sill.1 = P$c1, partial.sill.2 = P$c2,
            sill = P$c0 + P$c1 + P$c2, timescale.1 = P$phi1,
            timescale.2 = P$phi2)),
    power = list(
        coef = c("log_c0", "log_b", "logit_a"),
        g = function(h, P) P$c0 + P$b * h^P$alpha,
        report = function(P) c(nugget = P$c0, scale = P$b,
            exponent = P$alpha)),
    linear = list(
        coef = c("log_c0", "log_b"),
        g = function(h, P) P$c0 + P$b * h,
        report = function(P) c(nugget = P$c0, slope = P$b))
)

# Natural parameters from the coefficients
.natural_parameters <- function(b, form){
    names(b) <- .VARIOGRAM_FORMS[[form]]$coef
    P <- list()
    for( n in names(b) ){
        if( n == "logit_a" ){
            P$alpha <- 2 * stats::plogis(b[[n]])
        } else{
            P[[sub("^log_", "", n)]] <- exp(b[[n]])
        }
    }
    if( !is.null(P$r) ){
        P$phi2 <- P$phi1 * (1 + P$r)
    }
    return(P)
}

.variogram_form <- function(h, b, form){
    .VARIOGRAM_FORMS[[form]]$g(h, .natural_parameters(b, form))
}

# The longest lag reached by at least max(5, 20%) of the subjects
.max_lag <- function(pairs){
    hs <- sort(tapply(pairs$h, pairs$subject, max), decreasing = TRUE)
    k <- max(5L, ceiling(0.2 * length(hs)))
    if( length(hs) < k ){
        return(max(hs))
    }
    return(unname(hs[k]))
}

# Weights that make every subject count equally, with mean one
.subject_weights <- function(subject){
    w <- as.numeric(1 / table(subject)[subject])
    return(w / mean(w))
}

# Starting values from the lag-binned empirical variogram
.start_values <- function(pairs, form){
    p <- pairs[pairs$h > 0, ]
    n_bins <- min(10L, max(3L, floor(nrow(p) / 20)))
    breaks <- unique(stats::quantile(log(p$h),
        probs = seq(0, 1, length.out = n_bins + 1L)))
    breaks <- exp(breaks)
    breaks[1L] <- breaks[1L] * 0.999
    bin <- cut(p$h, breaks, include.lowest = TRUE)
    w <- as.numeric(1 / table(p$subject)[p$subject])
    h_bin <- tapply(w * p$h, bin, sum) / tapply(w, bin, sum)
    d_bin <- tapply(w * p$d, bin, sum) / tapply(w, bin, sum)
    d_bin <- d_bin[order(h_bin)]
    d_bin <- d_bin[!is.na(d_bin)]
    d_low <- d_bin[[1L]]
    sill <- mean(utils::tail(d_bin, min(3L, length(d_bin))))
    c0 <- max(0.5 * d_low, 1e-6 * sill)
    c1 <- max(sill - c0, 0.1 * sill)
    hmin <- max(min(p$h), 1e-3 * max(p$h))
    hmax <- max(p$h)
    phis <- exp(seq(log(hmin), log(hmax), length.out = 6L))
    hmed <- stats::median(p$h)
    starts <- switch(form,
        exponential = lapply(phis, function(ph){
            c(log(c0), log(c1), log(ph))
        }),
        two.timescale = unlist(lapply(phis[-length(phis)], function(ph){
            lapply(c(5, 30), function(r){
                c(log(c0), log(c1 / 2), log(ph), log(c1 / 2), log(r - 1))
            })
        }), recursive = FALSE),
        power = lapply(c(0.3, 1, 1.7), function(a){
            c(log(c0), log(c1 / hmed^a), stats::qlogis(a / 2))
        }),
        linear = list(c(log(c0), log(c1 / hmed))))
    return(starts)
}

# Weighted least squares with variance proportional to g(h)^2, by iterated
# weights from several starting values; the fit with the smallest deviance
# is kept
.fit_variogram <- function(pairs, form, starts = NULL, maxit = 30L){
    w <- .subject_weights(pairs$subject)
    d <- pairs$d
    h <- pairs$h
    floor_g <- 1e-8 * mean(d)
    if( is.null(starts) ){
        starts <- .start_values(pairs, form)
    }
    fits <- lapply(starts, function(b){
        g <- pmax(.variogram_form(h, b, form), floor_g)
        if( any(!is.finite(g)) ){
            return(list(beta = b, deviance = Inf, converged = FALSE))
        }
        converged <- FALSE
        for( it in seq_len(maxit) ){
            sw <- sqrt(w) / g
            res_fn <- function(bb){
                r <- sw * (d - .variogram_form(h, bb, form))
                r[!is.finite(r)] <- 1e10
                r
            }
            b_new <- .levenberg_marquardt(b, res_fn)
            g_new <- pmax(.variogram_form(h, b_new, form), floor_g)
            if( any(!is.finite(g_new)) ){
                break
            }
            change <- max(abs(b_new - b) / (abs(b) + 0.1))
            b <- b_new
            g <- g_new
            if( change < 1e-5 ){
                converged <- TRUE
                break
            }
        }
        dev <- sum(w * (d / g - 1)^2)
        list(beta = b, deviance = if( is.finite(dev) ) dev else Inf,
            converged = converged)
    })
    dev <- vapply(fits, function(f) f$deviance, numeric(1L))
    if( all(!is.finite(dev)) ){
        stop("The fit of the ", form, " form failed from all starting ",
            "values.", call. = FALSE)
    }
    best <- fits[[which.min(dev)]]
    res <- list(
        beta = best$beta, form = form, converged = best$converged,
        deviance = best$deviance, pairs = pairs, w = w,
        fitted = .variogram_form(h, best$beta, form))
    return(res)
}

# Levenberg-Marquardt minimisation of the sum of squared residuals
# (MINPACK), with quasi-Newton minimisation if it fails
.levenberg_marquardt <- function(b, res_fn){
    res <- tryCatch(minpack.lm::nls.lm(par = b, fn = res_fn,
        control = minpack.lm::nls.lm.control(maxiter = 200L)),
        error = function(e) NULL)
    if( !is.null(res) && all(is.finite(res$par)) ){
        return(res$par)
    }
    res <- stats::optim(b, function(bb) sum(res_fn(bb)^2), method = "BFGS",
        control = list(maxit = 500L))
    return(res$par)
}

# Central differences of f with respect to each element of b
.numerical_jacobian <- function(f, b, step = 1e-6){
    J <- vapply(seq_along(b), function(j){
        e <- step * max(1, abs(b[[j]]))
        bp <- b
        bm <- b
        bp[[j]] <- b[[j]] + e
        bm[[j]] <- b[[j]] - e
        (f(bp) - f(bm)) / (2 * e)
    }, numeric(length(f(b))))
    return(matrix(J, ncol = length(b)))
}

# Cluster-robust (CR1) covariance of the coefficients, subjects as clusters
.fit_vcov <- function(fit){
    p <- fit$pairs
    J <- .numerical_jacobian(
        function(b) .variogram_form(p$h, b, fit$form), fit$beta, 1e-5)
    wv <- fit$w / pmax(fit$fitted, 1e-8 * mean(p$d))^2
    A <- crossprod(J, J * wv)
    U <- rowsum(J * (wv * (p$d - fit$fitted)), p$subject)
    G <- nrow(U)
    Ai <- tryCatch(solve(A), error = function(e) .generalised_inverse(A))
    V <- Ai %*% crossprod(U) %*% Ai * G / (G - 1)
    return(V)
}

.generalised_inverse <- function(A){
    s <- svd(A)
    pos <- s$d > max(dim(A)) * max(s$d) * .Machine$double.eps
    s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
}

# Reported parameters with intervals on the log scale (delta method)
.variogram_parameters <- function(fit, conf.level){
    rep_fn <- .VARIOGRAM_FORMS[[fit$form]]$report
    est <- rep_fn(.natural_parameters(fit$beta, fit$form))
    V <- .fit_vcov(fit)
    J <- .numerical_jacobian(
        function(b) log(rep_fn(.natural_parameters(b, fit$form))), fit$beta)
    se <- sqrt(pmax(rowSums((J %*% V) * J), 0))
    G <- length(unique(fit$pairs$subject))
    q <- stats::qt(1 - (1 - conf.level) / 2, G - 1)
    h <- fit$pairs$h
    min_lag <- min(h[h > 0])
    max_lag <- max(h)
    lower <- est * exp(-q * se)
    upper <- est * exp(q * se)
    timescale <- grepl("^timescale", names(est))
    res <- data.frame(
        type = "parameter", form = fit$form, parameter = names(est),
        estimate = unname(est), lower = unname(lower), upper = unname(upper),
        se = se, df = G - 1, p.value = NA_real_,
        identified = ifelse(timescale,
            is.finite(se) & lower >= min_lag & upper <= max_lag, NA),
        n.subjects = G, n.pairs = nrow(fit$pairs), min.lag = min_lag,
        max.lag = max_lag)
    rownames(res) <- NULL
    return(res)
}

# Fits the form to each level of split.by over a common lag range, with
# percentile intervals from a bootstrap over subjects within each level, and
# compares each level with the first by the ratio of each parameter
.compare_variogram_models <- function(
        pairs, form, conf.level, bootstrap, BPPARAM){
    lev <- levels(droplevels(pairs$level))
    by_level <- split(pairs, pairs$level)[lev]
    max_lag <- min(vapply(by_level, .max_lag, numeric(1L)))
    by_level <- lapply(by_level, function(p) p[p$h <= max_lag, , drop = FALSE])
    cv <- NULL
    if( length(form) > 1L ){
        cv <- .cv_forms(do.call(rbind, unname(by_level)), form)
        form <- cv$form[cv$chosen]
    }
    a <- (1 - conf.level) / 2
    fits <- lapply(lev, function(g){
        p <- by_level[[g]]
        fit <- .fit_variogram(p, form)
        if( !fit$converged ){
            warning("The fit of the ", form, " form did not converge for '",
                g, "'.", call. = FALSE)
        }
        est <- .VARIOGRAM_FORMS[[form]]$report(
            .natural_parameters(fit$beta, form))
        boot <- .bootstrap_subjects(p, form, bootstrap, BPPARAM)
        if( any(is.na(boot[, 1L])) ){
            message(sum(is.na(boot[, 1L])), " of ", bootstrap,
                " resamples of '", g, "' failed and were left out.")
        }
        ci <- apply(boot, 2L, stats::quantile, probs = c(a, 1 - a),
            na.rm = TRUE)
        h <- p$h
        min_lag <- min(h[h > 0])
        timescale <- grepl("^timescale", names(est))
        res <- data.frame(
            level = g, type = "parameter", form = form,
            parameter = names(est), estimate = unname(est),
            lower = ci[1L, ], upper = ci[2L, ], se = NA_real_, df = NA_real_,
            p.value = NA_real_,
            identified = ifelse(timescale,
                ci[1L, ] >= min_lag & ci[2L, ] <= max_lag, NA),
            n.subjects = length(unique(p$subject)), n.pairs = nrow(p),
            min.lag = min_lag, max.lag = max_lag)
        list(table = res, boot = boot)
    })
    comp <- lapply(seq_along(lev)[-1L], function(i){
        e1 <- fits[[1L]]$table
        e2 <- fits[[i]]$table
        ok <- stats::complete.cases(fits[[1L]]$boot) &
            stats::complete.cases(fits[[i]]$boot)
        lr <- log(fits[[i]]$boot[ok, , drop = FALSE] /
            fits[[1L]]$boot[ok, , drop = FALSE])
        # A parameter estimated at 0 in both resamples has no ratio
        lr[!is.finite(lr)] <- NA
        p <- apply(lr, 2L, function(v){
            min(1, 2 * min(mean(v <= 0, na.rm = TRUE),
                mean(v >= 0, na.rm = TRUE)))
        })
        # Timescales are compared only when identified in both levels
        p[!is.na(e1$identified) & !(e1$identified & e2$identified)] <- NA
        data.frame(
            level = paste(lev[[i]], "/", lev[[1L]]), type = "ratio",
            form = form, parameter = e1$parameter,
            estimate = e2$estimate / e1$estimate,
            lower = exp(apply(lr, 2L, stats::quantile, probs = a,
                na.rm = TRUE)),
            upper = exp(apply(lr, 2L, stats::quantile, probs = 1 - a,
                na.rm = TRUE)),
            se = NA_real_, df = NA_real_, p.value = unname(p),
            identified = NA, n.subjects = e1$n.subjects + e2$n.subjects,
            n.pairs = e1$n.pairs + e2$n.pairs,
            min.lag = max(e1$min.lag, e2$min.lag), max.lag = max_lag)
    })
    res <- do.call(rbind, c(lapply(fits, function(f) f$table), comp))
    rownames(res) <- NULL
    attr(res, "cv") <- cv
    return(res)
}

# Parameters refitted to resamples of the subjects, with replacement; a
# subject drawn more than once enters as separate subjects. NA when the fit
# fails or does not converge.
.bootstrap_subjects <- function(pairs, form, bootstrap, BPPARAM){
    by_subject <- split(seq_len(nrow(pairs)), pairs$subject)
    n_par <- length(.VARIOGRAM_FORMS[[form]]$report(.natural_parameters(
        rep(0, length(.VARIOGRAM_FORMS[[form]]$coef)), form)))
    res <- BiocParallel::bplapply(seq_len(bootstrap), function(b){
        draw <- sample(names(by_subject), length(by_subject), replace = TRUE)
        rows <- by_subject[draw]
        pb <- pairs[unlist(rows), , drop = FALSE]
        pb$subject <- rep(paste0(draw, "#", seq_along(draw)), lengths(rows))
        fit <- tryCatch(.fit_variogram(pb, form), error = function(e) NULL)
        if( is.null(fit) || !fit$converged ){
            return(rep(NA_real_, n_par))
        }
        .VARIOGRAM_FORMS[[form]]$report(.natural_parameters(fit$beta, form))
    }, BPPARAM = BPPARAM)
    res <- do.call(rbind, res)
    return(res)
}

# Ten-fold cross-validation over subjects. Single pairs are too noisy to
# separate forms, so the held-out subjects' variogram is binned at quantiles
# of the lags and compared with the fitted variogram by the mean squared
# relative difference. The simplest form within one standard error of the
# best is chosen.
.cv_forms <- function(pairs, forms, K = 10L, n_bins = 12L){
    subj <- unique(pairs$subject)
    K <- min(K, length(subj))
    fold <- stats::setNames(sample(rep_len(seq_len(K), length(subj))), subj)
    hp <- pairs$h[pairs$h > 0]
    breaks <- unique(stats::quantile(hp, seq(0, 1, length.out = n_bins + 1L)))
    breaks[1L] <- 0
    losses <- vapply(forms, function(form){
        vapply(seq_len(K), function(k){
            test_s <- names(fold)[fold == k]
            train <- pairs[!pairs$subject %in% test_s, , drop = FALSE]
            test <- pairs[pairs$subject %in% test_s, , drop = FALSE]
            f <- tryCatch(.fit_variogram(train, form),
                error = function(e) NULL)
            if( is.null(f) ){
                return(NA_real_)
            }
            g <- .variogram_form(test$h, f$beta, form)
            bin <- as.integer(cut(test$h, breaks, include.lowest = TRUE))
            w <- as.numeric(1 / table(test$subject)[test$subject])
            dbar <- tapply(w * test$d, bin, sum) / tapply(w, bin, sum)
            gbar <- tapply(w * g, bin, sum) / tapply(w, bin, sum)
            mean((dbar / gbar - 1)^2)
        }, numeric(1L))
    }, numeric(K))
    losses <- matrix(losses, ncol = length(forms),
        dimnames = list(NULL, forms))
    m <- colMeans(losses, na.rm = TRUE)
    n_coef <- vapply(forms, function(f) length(.VARIOGRAM_FORMS[[f]]$coef),
        integer(1L))
    if( all(is.na(m)) ){
        stop("No form could be fitted in the cross-validation.",
            call. = FALSE)
    }
    best <- which.min(m)
    diff <- losses - losses[, best]
    se_diff <- apply(diff, 2L, stats::sd, na.rm = TRUE) /
        sqrt(colSums(!is.na(diff)))
    se_diff[best] <- 0
    eligible <- names(m)[m - m[best] <= se_diff & !is.na(m)]
    eligible <- eligible[n_coef[eligible] == min(n_coef[eligible])]
    choice <- eligible[which.min(m[eligible])]
    res <- data.frame(
        form = forms, loss = unname(m), difference = unname(m - m[best]),
        se = unname(se_diff), n.coefficients = unname(n_coef),
        chosen = forms == choice)
    return(res)
}
