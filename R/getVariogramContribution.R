#' @name
#' getVariogramContribution
#'
#' @export
#'
#' @title
#' Feature contributions to community change at stated time lags
#'
#' @description
#' Splits \eqn{D(L)}, the expected dissimilarity between two samples of the
#' same subject taken \eqn{L} time units apart, into the parts carried by
#' each feature, for the squared Euclidean distance. This shows which
#' features carry the change at a lag and, between groups, which carry a
#' difference in change.
#'
#' @details
#' The squared Euclidean distance between two samples is the sum over
#' features of their squared differences. \eqn{D(L)} is estimated as in
#' \code{\link[miaTime:getVariogram]{getVariogram()}}, by a local regression
#' whose estimate is linear in the dissimilarities, so it splits exactly into
#' per-feature parts: the contribution of a feature is the same estimate
#' computed from that feature's squared differences, and the contributions
#' sum to \eqn{D(L)}. The window, the weights (every subject counts equally),
#' the support rules and \code{time.point} are those of
#' \code{getVariogram()}. Transform the abundances first, for example with
#' \code{"hellinger"} or \code{"clr"} in
#' \code{\link[mia:transformAssay]{transformAssay()}}, so that the distance
#' is the Hellinger or Aitchison distance. The local regression can give a
#' small negative contribution to a feature that hardly changes.
#'
#' Standard errors are cluster-robust by subject (CR1), with \eqn{G - 1}
#' degrees of freedom, and intervals are on the scale of the contribution.
#'
#' With \code{split.by}, contributions are estimated for each level, and each
#' level is compared with the first by the difference of the contributions.
#' For a subject-level variable the levels are independent and the
#' difference has a Welch–Satterthwaite t test. For a sample-level state,
#' each pair takes the level of its earlier sample, and the difference is the
#' coefficient of the level indicator in one regression on the pairs of both
#' levels, as for the ratio in \code{getVariogram()}. P-values are not
#' adjusted for multiple comparisons.
#'
#' @return
#' A \code{DataFrame} with one row per lag, feature (and level of
#' \code{split.by}, and comparison) and the columns
#' \describe{
#'   \item{lag}{The lag \eqn{L}.}
#'   \item{\code{split.by}}{When \code{split.by} is given: the level, or the
#'   comparison (\code{"B - A"}).}
#'   \item{type}{\code{"contribution"}, or \code{"difference"} for
#'   comparisons between levels.}
#'   \item{feature}{The feature.}
#'   \item{estimate, lower, upper}{The contribution to \eqn{D(L)}, or the
#'   difference of contributions, with the confidence interval.}
#'   \item{share}{For contributions, the share of \eqn{D(L)}; the shares of a
#'   lag sum to one.}
#'   \item{se}{Standard error.}
#'   \item{df}{Degrees of freedom of the t distribution.}
#'   \item{p.value}{For comparisons, the p-value of the test that the
#'   difference is 0.}
#'   \item{n.subjects, n.pairs}{Subjects and pairs in the lag window.}
#'   \item{supported}{Whether the design supports \eqn{D(L)} at this lag.}
#'   \item{time.point}{The time at which \eqn{D(L)} was estimated, when
#'   \code{time.point} is given.}
#' }
#'
#' @inheritParams getVariogram
#'
#' @param time.lag \code{NULL} or \code{numeric vector}. The lags \eqn{L}, in
#' the units of \code{time.col}. If \code{NULL}, ten lags between the
#' shortest and longest supported lag. (Default: \code{NULL})
#'
#' @param ... additional arguments.
#' \itemize{
#'   \item \code{altexp}: \code{Character scalar} or \code{integer scalar}.
#'   Specifies an alternative experiment containing the input data.
#' }
#'
#' @examples
#' library(miaTime)
#'
#' data(hitchip1006)
#' tse <- transformAssay(hitchip1006, method = "hellinger")
#'
#' # Features carrying change over 1 time unit, as squared Hellinger
#' # distance
#' res <- getVariogramContribution(
#'     tse,
#'     time.col = "time",
#'     group = "subject",
#'     assay.type = "hellinger",
#'     time.lag = 1
#'     )
#' res <- res[order(-res$share), ]
#' head(res)
#'
#' @seealso
#' \code{\link[miaTime:getVariogram]{getVariogram()}}
#'
NULL

#' @rdname getVariogramContribution
#' @export
setMethod("getVariogramContribution",
    signature = c(x = "SummarizedExperiment"),
    function(
        x,
        time.col,
        group,
        assay.type = "counts",
        time.lag = NULL,
        split.by = NULL,
        time.point = NULL,
        window = 2,
        conf.level = 0.95,
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
        temp <- .check_input(time.lag, list(NULL, "numeric vector"))
        if( !is.null(time.lag) && (anyNA(time.lag) || any(time.lag <= 0)) ){
            stop("'time.lag' must be positive.", call. = FALSE)
        }
        temp <- .check_input(
            split.by, list(NULL, "character scalar"), colnames(colData(x)))
        .check_split_by_name(
            split.by, c(.VARIOGRAM_COLUMNS, "feature", "share"))
        if( !(is.null(time.point) || identical(time.point, "median") ||
                (is.numeric(time.point) && !anyNA(time.point))) ){
            stop("'time.point' must be NULL, numeric or \"median\".",
                call. = FALSE)
        }
        temp <- .check_input(window, list("numeric scalar"))
        if( window <= 1 ){
            stop("'window' must be greater than 1.", call. = FALSE)
        }
        temp <- .check_input(conf.level, list("numeric scalar"))
        if( conf.level <= 0 || conf.level >= 1 ){
            stop("'conf.level' must be between 0 and 1.", call. = FALSE)
        }
        ########################### INPUT CHECK END ############################
        pairs <- .get_within_subject_pairs(
            x, time.col, group, assay.type, "euclidean", split.by, TRUE,
            FALSE, list())
        if( is.null(time.lag) ){
            time.lag <- .default_lags(pairs)
        }
        mat <- t(assay(x, assay.type))[attr(pairs, "samples"), , drop = FALSE]
        feat <- rownames(x)
        if( is.null(feat) ){
            feat <- paste0("feature", seq_len(nrow(x)))
        }
        colnames(mat) <- feat
        res <- .variogram_contribution(
            pairs, mat, time.lag, split.by, time.point, window, conf.level)
        res <- DataFrame(res, check.names = FALSE)
        return(res)
    }
)

################################ HELP FUNCTIONS ################################

# Contributions for all lags and levels of split.by, and differences from the
# first level
.variogram_contribution <- function(
        pairs, mat, lags, split.by, time.point, window, conf.level){
    a0 <- .get_at_time(pairs, lags, time.point, window)
    if( is.null(split.by) ){
        res <- .contribution(pairs, mat, lags, window, a0)
        res <- .format_contribution(res, time.point, conf.level)
        return(res)
    }
    lev <- levels(droplevels(pairs$level))
    est <- lapply(lev, function(g){
        .contribution(pairs[pairs$level == g, ], mat, lags, window, a0)
    })
    comp <- lapply(seq_along(lev)[-1L], function(i){
        if( isTRUE(attr(pairs, "pair.level")) ){
            .contribution_pair_levels(
                pairs, mat, lev[[1L]], lev[[i]], lags, window, a0,
                est[[1L]], est[[i]])
        } else{
            .contribution_difference(est[[1L]], est[[i]])
        }
    })
    est <- lapply(seq_along(lev), function(i){
        cbind(level = lev[[i]],
            .format_contribution(est[[i]], time.point, conf.level))
    })
    comp <- lapply(seq_along(comp), function(i){
        lab <- paste(lev[[i + 1L]], "-", lev[[1L]])
        cbind(level = lab,
            .format_contribution(comp[[i]], time.point, conf.level))
    })
    res <- do.call(rbind, c(est, comp))
    res <- res[order(res$lag), , drop = FALSE]
    rownames(res) <- NULL
    names(res)[names(res) == "level"] <- split.by
    res <- res[, c("lag", split.by, setdiff(names(res), c("lag", split.by)))]
    return(res)
}

# Squared differences of every feature for the pairs in rows
.squared_differences <- function(pairs, mat, rows){
    (mat[pairs$k[rows], , drop = FALSE] -
        mat[pairs$l[rows], , drop = FALSE])^2
}

# Each feature's part of the local estimate of D(L), with its CR1 standard
# error: the estimate is a weighted sum of the dissimilarities, so it is the
# same weighted sum of each feature's squared differences
.contribution <- function(pairs, mat, lags, window, a0){
    k <- max(5L, ceiling(0.2 * length(unique(pairs$subject))))
    res <- lapply(seq_along(lags), function(i){
        f <- .local_design(pairs, lags[[i]], a0[[i]], k, window)
        out <- data.frame(
            lag = lags[[i]], feature = colnames(mat), estimate = NA_real_,
            se = NA_real_, df = NA_real_, n.subjects = f$G,
            n.pairs = length(f$rows), supported = FALSE,
            time.point = a0[[i]])
        if( !f$ok ){
            return(out)
        }
        out[, c("estimate", "se")] <- .coefficient_by_feature(
            f, .squared_differences(pairs, mat, f$rows), 1L)
        out$df <- f$G - 1
        out$supported <- sum(out$estimate) > 0
        return(out)
    })
    res <- do.call(rbind, res)
    return(res)
}

# Coefficient j of the local regression of every column of Y, with its CR1
# standard error
.coefficient_by_feature <- function(f, Y, j){
    a <- f$Ai[j, ]
    cw <- drop(f$X %*% a) * f$w
    B <- f$Ai %*% crossprod(f$X, f$w * Y)
    R <- Y - f$X %*% B
    S <- rowsum(cw * R, f$subj)
    G <- f$G
    data.frame(estimate = B[j, ], se = sqrt(colSums(S^2) * G / (G - 1)))
}

# Difference of contributions between two independent levels, with a
# Welch-Satterthwaite t test
.contribution_difference <- function(e1, e2){
    est <- e2$estimate - e1$estimate
    se <- sqrt(e1$se^2 + e2$se^2)
    df <- se^4 / (e1$se^4 / e1$df + e2$se^4 / e2$df)
    ok <- e1$supported & e2$supported
    data.frame(
        lag = e1$lag, feature = e1$feature, estimate = est, se = se, df = df,
        n.subjects = e1$n.subjects + e2$n.subjects,
        n.pairs = e1$n.pairs + e2$n.pairs, supported = ok,
        time.point = e1$time.point,
        p.value = ifelse(ok, 2 * stats::pt(-abs(est / se), df), NA_real_),
        type = "difference")
}

# Difference of contributions between two levels of a sample-level state:
# the coefficient of the level indicator in one regression on the pairs of
# both levels
.contribution_pair_levels <- function(
        pairs, mat, lev1, lev2, lags, window, a0, e1, e2){
    p <- pairs[pairs$level %in% c(lev1, lev2), , drop = FALSE]
    res <- lapply(seq_along(lags), function(i){
        f <- .pair_level_design(p, lev2, lags[[i]], a0[[i]], window)
        out <- data.frame(
            lag = lags[[i]], feature = colnames(mat), estimate = NA_real_,
            se = NA_real_, df = NA_real_, n.subjects = f$G,
            n.pairs = length(f$rows), supported = FALSE,
            time.point = a0[[i]], p.value = NA_real_, type = "difference")
        ok <- f$ok && all(e1$supported[e1$lag == lags[[i]]]) &&
            all(e2$supported[e2$lag == lags[[i]]])
        if( !ok ){
            return(out)
        }
        out[, c("estimate", "se")] <- .coefficient_by_feature(
            f, .squared_differences(p, mat, f$rows), 2L)
        out$df <- f$G - 1
        out$supported <- TRUE
        out$p.value <- 2 * stats::pt(-abs(out$estimate / out$se), out$df)
        return(out)
    })
    res <- do.call(rbind, res)
    return(res)
}

# Intervals, shares, and columns in a fixed order
.format_contribution <- function(res, time.point, conf.level){
    if( is.null(res$type) ){
        res$type <- "contribution"
        res$p.value <- NA_real_
        res$share <- res$estimate / stats::ave(res$estimate, res$lag, FUN = sum)
    } else{
        res$share <- NA_real_
    }
    q <- stats::qt(1 - (1 - conf.level) / 2, res$df)
    res$lower <- res$estimate - q * res$se
    res$upper <- res$estimate + q * res$se
    res[!res$supported, c("estimate", "lower", "upper", "se", "share",
        "p.value")] <- NA_real_
    cols <- c("lag", "type", "feature", "estimate", "lower", "upper", "share",
        "se", "df", "p.value", "n.subjects", "n.pairs", "supported")
    if( !is.null(time.point) ){
        cols <- c(cols, "time.point")
    }
    res <- res[, cols]
    return(res)
}
