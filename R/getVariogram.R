#' @name
#' getVariogram
#'
#' @export
#'
#' @title
#' Community change at stated time lags (temporal variogram)
#'
#' @description
#' Estimates \eqn{D(L)}, the expected dissimilarity between two samples of the
#' same subject taken \eqn{L} time units apart, from all within-subject sample
#' pairs. Unlike divergence to the previous sample, \eqn{D(L)} refers to a
#' stated lag, so it can be compared between subjects, groups and studies
#' whose sampling schedules differ.
#'
#' @details
#' Every pair of samples from the same subject measures how much the
#' community changed over that pair's time lag \eqn{h}. The temporal
#' variogram \eqn{g(h)} is the expected dissimilarity as a function of the
#' lag, and \eqn{D(L) = g(L)}. Divergence to the previous sample measures
#' change at whatever intervals the schedule happens to have, so it differs
#' between designs even when the dynamics are the same.
#'
#' \eqn{D(L)} is estimated locally. Pairs with lags in
#' \eqn{[L / window, L \times window]} are used, and the dissimilarities are
#' regressed on \eqn{u = \log(h / L)} and \eqn{u^2} by weighted least
#' squares; the intercept estimates \eqn{D(L)}. The pairs of each subject are
#' weighted so that every subject counts equally. Pairs of the same subject
#' share samples and are correlated, so the standard error is cluster-robust
#' by subject (CR1), and the interval for \eqn{\log D(L)} uses a t
#' distribution with \eqn{G - 1} degrees of freedom, where \eqn{G} is the
#' number of subjects with pairs in the window.
#'
#' \eqn{D(L)} is reported only at lags the design supports: at least
#' \eqn{\max(5, 20\%)} of the subjects with two or more samples must have
#' pairs in the window, and at least three subjects must have pairs at lags
#' up to \eqn{L} and three at lags from \eqn{L}, so that \eqn{D(L)} is
#' interpolated, not extrapolated. Elsewhere the estimate is \code{NA} and
#' \code{supported} is \code{FALSE}.
#'
#' When \code{time.lag = NULL}, ten lags are spread evenly on the log scale
#' between the shortest and longest supported lag of the whole data: the
#' smallest and largest lag reached by at least \eqn{\max(5, 20\%)} of the
#' subjects.
#'
#' \strong{Groups.} With \code{split.by}, \eqn{D(L)} is estimated separately
#' for each level of a subject-level variable, and each level is compared
#' with the first level by the ratio of \eqn{D(L)}, with a
#' Welch–Satterthwaite t test of the log ratio. The comparison is made only
#' at lags supported in both groups. P-values are not adjusted for multiple
#' comparisons.
#'
#' \strong{Developing communities.} When change at a given lag depends on
#' the time since a common origin (age in infant cohorts, time since an
#' intervention), groups whose samples were taken at different ages differ
#' in \eqn{D(L)} even with the same development. \code{time.point} estimates
#' \eqn{D(L)} at a stated time \eqn{a}: only pairs whose midpoint \eqn{m}
#' satisfies \eqn{(1 + m) / (1 + a)} in \eqn{[1/2, 2]} are used, the
#' regression adds \eqn{z = \log(1 + m) - \log(1 + a)} and \eqn{z^2}, and
#' at least three subjects must have midpoints up to \eqn{a} and three from
#' \eqn{a}. With \code{time.point = "median"}, \eqn{a} is the median midpoint of
#' all pairs in the lag window, pooled over the levels of \code{split.by}, so
#' that every group is evaluated at the same time. Time must then be
#' measured from the common origin and be non-negative.
#'
#' \strong{Dissimilarity.} Any dissimilarity of
#' \code{\link[mia:getDissimilarity]{getDissimilarity()}} can be used. The
#' squared Euclidean distance between Hellinger-transformed abundances
#' (\code{method = "euclidean"}, \code{squared = TRUE} on a
#' \code{"hellinger"} assay) adds technical noise to the dissimilarity
#' additively, so the noise does not change the shape of the variogram, only
#' its level. Dissimilarities depend on sequencing depth; rarefy the counts
#' first (\code{\link[mia:rarefyAssay]{rarefyAssay()}}) when depth differs
#' between samples or groups.
#'
#' @references
#' Diggle P.J. (1988) An approach to the analysis of repeated measurements.
#' Biometrics 44, 959–971. doi: 10.2307/2531727
#'
#' Lloyd-Price J., et al. (2019) Multi-omics of the gut microbial ecosystem in
#' inflammatory bowel diseases. Nature 569, 655–662.
#' doi: 10.1038/s41586-019-1237-9
#'
#' Björk J.R., et al. (2022) Synchrony and idiosyncrasy in the gut microbiome
#' of wild baboons. Nature Ecology & Evolution 6, 955–964.
#' doi: 10.1038/s41559-022-01773-4
#'
#' @return
#' \code{getVariogram} returns a \code{DataFrame} with one row per
#' lag (and level of \code{split.by}, and comparison) and the columns
#' \describe{
#'   \item{lag}{The lag \eqn{L}.}
#'   \item{\code{split.by}}{When \code{split.by} is given: the level, or the
#'   comparison (\code{"B / A"}).}
#'   \item{type}{\code{"dissimilarity"} for estimates of \eqn{D(L)},
#'   \code{"ratio"} for comparisons between levels.}
#'   \item{estimate, lower, upper}{\eqn{D(L)} or the ratio of \eqn{D(L)},
#'   with the confidence interval.}
#'   \item{se}{Standard error on the log scale.}
#'   \item{df}{Degrees of freedom of the t distribution.}
#'   \item{p.value}{For comparisons, the p-value of the test that the ratio
#'   is 1.}
#'   \item{n.subjects, n.pairs}{Subjects and pairs in the lag window.}
#'   \item{supported}{Whether the design supports \eqn{D(L)} at this lag.}
#'   \item{time.point}{The time at which \eqn{D(L)} was estimated, when
#'   \code{time.point} is given.}
#' }
#' \code{addVariogram} returns \code{x} with this table in its
#' \code{metadata}.
#'
#' @inheritParams addBaselineDivergence
#'
#' @param group \code{Character scalar}. Specifies a name of the column from
#' \code{colData} that identifies the subject (or other unit followed over
#' time) of each sample. Pairs are formed within subjects.
#'
#' @param time.lag \code{NULL} or \code{numeric vector}. The lags \eqn{L}, in
#' the units of \code{time.col}, at which \eqn{D(L)} is estimated. If
#' \code{NULL}, ten lags between the shortest and longest supported lag.
#' (Default: \code{NULL})
#'
#' @param split.by \code{NULL} or \code{Character scalar}. Specifies a name
#' of a column from \code{colData} with one value per subject. \eqn{D(L)} is
#' estimated for each of its levels, and each level is compared with the
#' first. (Default: \code{NULL})
#'
#' @param time.point \code{NULL}, \code{numeric vector} or \code{"median"}. The
#' time at which \eqn{D(L)} is estimated, for communities whose change
#' depends on time since a common origin; one value, or one per lag. See
#' details. (Default: \code{NULL})
#'
#' @param window \code{Numeric scalar}. Pairs with lags within a factor of
#' \code{window} of \eqn{L} are used. (Default: \code{2})
#'
#' @param squared \code{Logical scalar}. Whether the dissimilarities are
#' squared. (Default: \code{FALSE})
#'
#' @param conf.level \code{Numeric scalar}. Confidence level of the
#' intervals. (Default: \code{0.95})
#'
#' @param name \code{Character scalar}. Specifies a name for storing the
#' results in \code{metadata}. (Default: \code{"variogram"})
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
#' # Change over 1 to 9 time units, as squared Hellinger distance
#' res <- getVariogram(
#'     tse,
#'     time.col = "time",
#'     group = "subject",
#'     assay.type = "hellinger",
#'     method = "euclidean",
#'     squared = TRUE,
#'     time.lag = c(1, 3, 6, 9)
#'     )
#' res
#'
#' # Compare the sexes
#' res <- getVariogram(
#'     tse,
#'     time.col = "time",
#'     group = "subject",
#'     assay.type = "hellinger",
#'     method = "euclidean",
#'     squared = TRUE,
#'     time.lag = c(1, 3),
#'     split.by = "sex"
#'     )
#' res[res$type == "ratio", ]
#'
#' # Store the results in metadata
#' tse <- addVariogram(
#'     tse,
#'     time.col = "time",
#'     group = "subject",
#'     assay.type = "hellinger",
#'     method = "euclidean",
#'     squared = TRUE
#'     )
#' metadata(tse)[["variogram"]]
#'
#' @seealso
#' \code{\link[miaTime:getStepwiseDivergence]{getStepwiseDivergence()}}
#'
NULL

#' @rdname getVariogram
#' @export
setMethod("addVariogram", signature = c(x = "SummarizedExperiment"),
    function(x, name = "variogram", ...){
        temp <- .check_input(name, list("character scalar"))
        x <- .check_and_get_altExp(x, ...)
        args <- c(list(x = x), list(...)[!names(list(...)) %in% c("altexp")])
        res <- do.call(getVariogram, args)
        x <- .add_values_to_metadata(x, name, res, ...)
        return(x)
    }
)

#' @rdname getVariogram
#' @export
setMethod("getVariogram", signature = c(x = "SummarizedExperiment"),
    function(
        x,
        time.col,
        group,
        assay.type = "counts",
        method = "bray",
        time.lag = NULL,
        split.by = NULL,
        time.point = NULL,
        window = 2,
        squared = FALSE,
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
        temp <- .check_input(method, list("character scalar"))
        temp <- .check_input(time.lag, list(NULL, "numeric vector"))
        if( !is.null(time.lag) && (anyNA(time.lag) || any(time.lag <= 0)) ){
            stop("'time.lag' must be positive.", call. = FALSE)
        }
        temp <- .check_input(
            split.by, list(NULL, "character scalar"), colnames(colData(x)))
        if( !(is.null(time.point) || identical(time.point, "median") ||
                (is.numeric(time.point) && !anyNA(time.point))) ){
            stop("'time.point' must be NULL, numeric or \"median\".",
                call. = FALSE)
        }
        if( is.numeric(time.point) &&
                !length(time.point) %in% c(1L, length(time.lag)) ){
            stop("'time.point' must have length 1 or the length of ",
                "'time.lag'.",
                call. = FALSE)
        }
        temp <- .check_input(window, list("numeric scalar"))
        if( window <= 1 ){
            stop("'window' must be greater than 1.", call. = FALSE)
        }
        temp <- .check_input(squared, list("logical scalar"))
        temp <- .check_input(conf.level, list("numeric scalar"))
        if( conf.level <= 0 || conf.level >= 1 ){
            stop("'conf.level' must be between 0 and 1.", call. = FALSE)
        }
        ########################### INPUT CHECK END ############################
        args <- list(...)[!names(list(...)) %in% c("altexp")]
        pairs <- .get_within_subject_pairs(
            x, time.col, group, assay.type, method, split.by, squared, args)
        if( is.null(time.lag) ){
            time.lag <- .default_lags(pairs)
        }
        res <- .variogram(
            pairs, time.lag, split.by, time.point, window, conf.level)
        res <- DataFrame(res, check.names = FALSE)
        return(res)
    }
)

################################ HELP FUNCTIONS ################################

# All pairs of samples of the same subject, with their lag, midpoint time and
# dissimilarity.
.get_within_subject_pairs <- function(
        x, time.col, group, assay.type, method, split.by, squared, args){
    time <- x[[time.col]]
    subject <- x[[group]]
    keep <- !is.na(time) & !is.na(subject)
    if( !is.null(split.by) ){
        keep <- keep & !is.na(x[[split.by]])
    }
    if( any(!keep) ){
        message(sum(!keep), " samples with missing values in '", time.col,
            "', '", group, "'",
            if( !is.null(split.by) ) paste0(" or '", split.by, "'"),
            " were removed.")
        x <- x[, keep]
        time <- time[keep]
        subject <- subject[keep]
    }
    subject <- as.character(subject)
    level <- NULL
    if( !is.null(split.by) ){
        level <- as.character(x[[split.by]])
        n_levels <- tapply(level, subject, function(v) length(unique(v)))
        if( any(n_levels > 1L) ){
            stop("'split.by' must have one value per subject ('", group,
                "').", call. = FALSE)
        }
        if( length(unique(level)) < 2L ){
            stop("'split.by' must have at least two levels.", call. = FALSE)
        }
        # The first level is the reference of the comparisons
        level <- factor(level, levels = levels(factor(x[[split.by]])))
    }
    mat <- t(assay(x, assay.type))
    idx_list <- split(seq_along(subject), subject)
    idx_list <- idx_list[lengths(idx_list) >= 2L]
    if( length(idx_list) == 0L ){
        stop("No subject has two or more samples.", call. = FALSE)
    }
    res <- lapply(idx_list, function(idx){
        d <- do.call(getDissimilarity, c(
            list(x = mat[idx, , drop = FALSE], method = method), args))
        d <- as.matrix(d)
        ut <- which(upper.tri(d), arr.ind = TRUE)
        k <- idx[ut[, 1L]]
        l <- idx[ut[, 2L]]
        data.frame(
            subject = subject[k], h = abs(time[l] - time[k]),
            mid = (time[k] + time[l]) / 2, d = d[ut])
    })
    pairs <- do.call(rbind, unname(res))
    if( squared ){
        pairs$d <- pairs$d^2
    }
    if( !is.null(split.by) ){
        pairs$level <- level[match(pairs$subject, subject)]
    }
    return(pairs)
}

# The k-th smallest of the subjects' shortest lags and the k-th largest of
# their longest lags, k = max(5, 20% of subjects).
.default_lags <- function(pairs){
    p <- pairs[pairs$h > 0, ]
    lo <- sort(tapply(p$h, p$subject, min))
    hi <- sort(tapply(p$h, p$subject, max), decreasing = TRUE)
    k <- max(5L, ceiling(0.2 * length(lo)))
    if( length(lo) < k ){
        stop("Too few subjects to estimate the variogram: pairs from at ",
            "least five subjects are needed.", call. = FALSE)
    }
    lo <- unname(lo[k])
    hi <- unname(hi[k])
    if( hi <= lo ){
        return(lo)
    }
    lags <- exp(seq(log(lo), log(hi), length.out = 10L))
    lags <- unique(signif(lags, 3L))
    return(lags)
}

# D(L) for all lags, levels of split.by, and comparisons with the first level
.variogram <- function(
        pairs, lags, split.by, time.point, window, conf.level){
    a0 <- .get_at_time(pairs, lags, time.point, window)
    if( is.null(split.by) ){
        res <- .local_variogram(pairs, lags, window, conf.level, a0)
        res <- .format_variogram(res, time.point, conf.level)
        return(res)
    }
    lev <- levels(droplevels(pairs$level))
    est <- lapply(lev, function(g){
        .local_variogram(
            pairs[pairs$level == g, ], lags, window, conf.level, a0)
    })
    comp <- lapply(est[-1L], function(e2){
        .compare_levels(est[[1L]], e2)
    })
    est <- lapply(seq_along(lev), function(i){
        cbind(level = lev[[i]],
            .format_variogram(est[[i]], time.point, conf.level))
    })
    comp <- lapply(seq_along(comp), function(i){
        lab <- paste(lev[[i + 1L]], "/", lev[[1L]])
        cbind(level = lab, .format_variogram(comp[[i]], time.point, conf.level))
    })
    res <- do.call(rbind, c(est, comp))
    res <- res[order(res$lag), , drop = FALSE]
    rownames(res) <- NULL
    names(res)[names(res) == "level"] <- split.by
    res <- res[, c("lag", split.by, setdiff(names(res), c("lag", split.by)))]
    return(res)
}

# The time at which D(L) is estimated: NA (no adjustment), the values given,
# or the median midpoint of all pairs in each lag window
.get_at_time <- function(pairs, lags, time.point, window){
    if( is.null(time.point) ){
        return(rep(NA_real_, length(lags)))
    }
    if( is.numeric(time.point) ){
        return(rep_len(as.numeric(time.point), length(lags)))
    }
    a0 <- vapply(lags, function(l){
        in_window <- pairs$h >= l / window & pairs$h <= l * window
        if( any(in_window) ) stats::median(pairs$mid[in_window]) else NA_real_
    }, numeric(1L))
    return(a0)
}

# Local quadratic estimate of D(L) with a cluster-robust (CR1) interval
.local_variogram <- function(pairs, lags, window, conf.level, a0){
    k <- max(5L, ceiling(0.2 * length(unique(pairs$subject))))
    res <- lapply(seq_along(lags), function(i){
        .local_variogram_at(pairs, lags[[i]], a0[[i]], k, window, conf.level)
    })
    res <- do.call(rbind, res)
    return(res)
}

.local_variogram_at <- function(pairs, l, a0, k, window, conf.level){
    rows <- which(pairs$h >= l / window & pairs$h <= l * window)
    if( !is.na(a0) ){
        # Local in time too: midpoints within a factor of two of the target
        # on the 1 + time scale
        r <- (1 + pmax(pairs$mid[rows], 0)) / (1 + max(a0, 0))
        rows <- rows[r >= 1 / 2 & r <= 2]
    }
    q <- pairs[rows, , drop = FALSE]
    subj <- q$subject
    ns <- table(subj)
    G <- length(ns)
    res <- data.frame(
        lag = l, estimate = NA_real_, se = NA_real_, df = NA_real_,
        n.subjects = G, n.pairs = nrow(q), supported = FALSE, time.point = a0)
    # Interpolation only: pairs from at least three subjects on each side
    below <- length(unique(subj[q$h <= l]))
    above <- length(unique(subj[q$h >= l]))
    if( G < max(k, 3L) || below < 3L || above < 3L ){
        return(res)
    }
    u <- log(q$h / l)
    X <- cbind(1, u, u^2)
    if( !is.na(a0) ){
        if( length(unique(subj[q$mid <= a0])) < 3L ||
                length(unique(subj[q$mid >= a0])) < 3L ){
            return(res)
        }
        z <- log1p(pmax(q$mid, 0)) - log1p(max(a0, 0))
        X <- cbind(X, z, z^2)
    }
    # Each subject's pairs in the window count as one subject
    w <- as.numeric(1 / ns[subj])
    A <- crossprod(X, X * w)
    if( qr(A)$rank < ncol(X) ){
        return(res)
    }
    Ai <- solve(A)
    b <- Ai %*% crossprod(X, w * q$d)
    r <- q$d - drop(X %*% b)
    U <- rowsum(X * (w * r), subj)
    V <- Ai %*% crossprod(U) %*% Ai * G / (G - 1)
    D <- b[[1L]]
    if( !is.finite(D) || D <= 0 ){
        return(res)
    }
    res$estimate <- D
    res$se <- sqrt(V[1L, 1L]) / D
    res$df <- G - 1
    res$supported <- TRUE
    return(res)
}

# Ratio of D(L) between two levels, with a Welch-Satterthwaite t test of the
# log ratio
.compare_levels <- function(e1, e2){
    lr <- log(e2$estimate) - log(e1$estimate)
    se <- sqrt(e1$se^2 + e2$se^2)
    df <- se^4 / (e1$se^4 / e1$df + e2$se^4 / e2$df)
    ok <- e1$supported & e2$supported
    data.frame(
        lag = e1$lag, estimate = exp(lr), se = se, df = df,
        n.subjects = e1$n.subjects + e2$n.subjects,
        n.pairs = e1$n.pairs + e2$n.pairs, supported = ok,
        time.point = e1$time.point,
        p.value = ifelse(ok, 2 * stats::pt(-abs(lr / se), df), NA_real_),
        type = "ratio")
}

# Confidence intervals on the log scale, and columns in a fixed order
.format_variogram <- function(res, time.point, conf.level){
    if( is.null(res$type) ){
        res$type <- "dissimilarity"
        res$p.value <- NA_real_
    }
    q <- stats::qt(1 - (1 - conf.level) / 2, res$df)
    res$lower <- res$estimate * exp(-q * res$se)
    res$upper <- res$estimate * exp(q * res$se)
    cols <- c("lag", "type", "estimate", "lower", "upper", "se", "df",
        "p.value", "n.subjects", "n.pairs", "supported")
    if( !is.null(time.point) ){
        cols <- c(cols, "time.point")
    }
    res <- res[, cols]
    return(res)
}
