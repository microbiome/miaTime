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
#' interpolated, not extrapolated. The window must also contain at least
#' twice as many subjects as regression coefficients (6, or 10 with
#' \code{time.point}), because the sandwich variance is singular when there
#' are as many subjects as coefficients. Elsewhere the estimate is \code{NA} and
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
#' With \code{permutations > 0}, each comparison also gets a permutation
#' p-value: the levels of \code{split.by} are reassigned at random among
#' the subjects of the two levels compared, keeping their numbers, and the
#' t statistic of the log ratio is recomputed. When the levels share a
#' sampling schedule, subjects are exchangeable under the null hypothesis
#' and the permutation test is exact, whereas the t test can reject slightly
#' too often in real data. When the schedules differ, the labels are not
#' exchangeable and the permutation p-value is not valid; use the t test.
#' Set a seed (\code{set.seed()}) for reproducible permutation p-values.
#'
#' \strong{Kinds of pairs.} When \code{split.by} varies within subjects (a
#' sample-level state, such as community type or the phase of an
#' intervention), each pair takes the level of its earlier sample, and
#' \eqn{D(L)} is estimated for the pairs of each level. A subject then
#' contributes pairs to several levels, so the levels are not independent.
#' Each level is compared with the first in one weighted regression on the
#' pairs of both levels in the lag window, with an indicator of the level,
#' the terms in \eqn{u} for each level and, with \code{time.point}, the
#' terms in \eqn{z} shared by the levels; the log ratio of \eqn{D(L)} is
#' \eqn{\log((b_0 + b_1) / b_0)}, where \eqn{b_0} is the intercept and
#' \eqn{b_1} the coefficient of the indicator. Its standard error is
#' cluster-robust by subject, by the delta method, with \eqn{G - 1} degrees
#' of freedom. Both levels must meet the support rules on their own. The
#' permutation test is not available, as levels are not assigned to
#' subjects.
#'
#' With \code{subset.by}, only pairs whose earlier sample is \code{TRUE} in
#' a logical column of \code{colData} are used, for example pairs starting
#' from a given community type. Subjects are then compared on those pairs,
#' and a sample-level state is compared within them.
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
#' Dissimilarities bounded by 1, such as Bray–Curtis and Jaccard, approach
#' their bound as communities share fewer features, so \eqn{D(L)} flattens
#' near 1 whatever the dynamics. With \code{bounded = TRUE} they are analysed
#' on the scale \eqn{-\log(1 - d)}, which is unbounded and close to \eqn{d}
#' for small dissimilarities. Pairs without shared features (\eqn{d = 1}) are
#' set to \eqn{d = 1 - 10^{-9}}.
#'
#' \strong{Relative change.} \eqn{D(L)} is on the scale of the
#' dissimilarity, so it cannot be compared between dissimilarities, studies
#' or body sites. With \code{relative = TRUE}, \eqn{D(L)} is also reported
#' relative to \eqn{B}, the mean dissimilarity between samples of different
#' subjects, where every pair of subjects counts equally: a value of 0.5
#' means that over \eqn{L} a subject's community changes half as much as it
#' differs from another subject's. \eqn{B} is computed from the samples of
#' all subjects, so it is the same for every level of \code{split.by}, and
#' it is treated as fixed.
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
#'   \item{p.value.permutation}{When \code{permutations > 0}: for
#'   comparisons, the permutation p-value (see details).}
#'   \item{n.subjects, n.pairs}{Subjects and pairs in the lag window.}
#'   \item{supported}{Whether the design supports \eqn{D(L)} at this lag.}
#'   \item{time.point}{The time at which \eqn{D(L)} was estimated, when
#'   \code{time.point} is given.}
#'   \item{relative}{When \code{relative = TRUE}: \eqn{D(L) / B}, for
#'   estimates of \eqn{D(L)}.}
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
#' of a column from \code{colData}, with one value per subject or a
#' sample-level state. \eqn{D(L)} is estimated for each of its levels, and
#' each level is compared with the first. (Default: \code{NULL})
#'
#' @param subset.by \code{NULL} or \code{Character scalar}. Specifies a name
#' of a logical column of \code{colData}; only pairs whose earlier sample is
#' \code{TRUE} are used. See details. (Default: \code{NULL})
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
#' @param bounded \code{Logical scalar}. Whether the dissimilarities are
#' bounded by 1 and analysed on the scale \eqn{-\log(1 - d)}. See details.
#' (Default: \code{FALSE})
#'
#' @param relative \code{Logical scalar}. Whether \eqn{D(L)} is also
#' reported relative to the dissimilarity between subjects. See details.
#' (Default: \code{FALSE})
#'
#' @param conf.level \code{Numeric scalar}. Confidence level of the
#' intervals. (Default: \code{0.95})
#'
#' @param permutations \code{Integer scalar}. The number of permutations
#' of \code{split.by} among subjects for a permutation p-value of each
#' comparison; \code{0} for none. See details. (Default: \code{0})
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
        subset.by = NULL,
        time.point = NULL,
        window = 2,
        squared = FALSE,
        bounded = FALSE,
        relative = FALSE,
        conf.level = 0.95,
        permutations = 0,
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
        .check_split_by_name(split.by, .VARIOGRAM_COLUMNS)
        temp <- .check_input(
            subset.by, list(NULL, "character scalar"), colnames(colData(x)))
        if( !is.null(subset.by) && !is.logical(x[[subset.by]]) ){
            stop("'subset.by' must specify a logical column of colData(x).",
                call. = FALSE)
        }
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
        temp <- .check_input(bounded, list("logical scalar"))
        temp <- .check_input(relative, list("logical scalar"))
        temp <- .check_input(conf.level, list("numeric scalar"))
        if( conf.level <= 0 || conf.level >= 1 ){
            stop("'conf.level' must be between 0 and 1.", call. = FALSE)
        }
        temp <- .check_input(permutations, list("numeric scalar"))
        if( permutations < 0 || permutations != round(permutations) ){
            stop("'permutations' must be a non-negative integer.",
                call. = FALSE)
        }
        ########################### INPUT CHECK END ############################
        args <- list(...)[!names(list(...)) %in% c("altexp")]
        pairs <- .get_within_subject_pairs(
            x, time.col, group, assay.type, method, split.by, squared,
            bounded, args)
        if( !is.null(subset.by) ){
            pairs <- .subset_pairs(pairs, x, time.col, subset.by)
        }
        if( is.null(time.lag) ){
            time.lag <- .default_lags(pairs)
        }
        if( isTRUE(attr(pairs, "pair.level")) && permutations > 0 ){
            stop("'permutations' needs 'split.by' with one value per ",
                "subject.", call. = FALSE)
        }
        res <- .variogram(
            pairs, time.lag, split.by, time.point, window, conf.level,
            permutations)
        if( relative ){
            B <- .between_subject_level(
                x, time.col, group, assay.type, method, split.by, squared,
                bounded, args)
            res$relative <- ifelse(
                res$type == "dissimilarity", res$estimate / B, NA_real_)
        }
        res <- DataFrame(res, check.names = FALSE)
        return(res)
    }
)

################################ HELP FUNCTIONS ################################

# Columns of the output; split.by becomes a column too, so it must not take
# one of these names
.VARIOGRAM_COLUMNS <- c(
    "lag", "type", "estimate", "lower", "upper", "se", "df", "p.value",
    "p.value.permutation", "n.subjects", "n.pairs", "supported",
    "time.point", "relative")

.check_split_by_name <- function(split.by, columns){
    if( !is.null(split.by) && split.by %in% columns ){
        stop("'split.by' must not be named '", split.by, "', the name of ",
            "an output column; rename the column of colData.", call. = FALSE)
    }
}

# All pairs of samples of the same subject, with their lag, midpoint time and
# dissimilarity.
.get_within_subject_pairs <- function(
        x, time.col, group, assay.type, method, split.by, squared, bounded,
        args){
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
    samples <- which(keep)
    subject <- as.character(subject)
    level <- NULL
    if( !is.null(split.by) ){
        level <- as.character(x[[split.by]])
        n_levels <- tapply(level, subject, function(v) length(unique(v)))
        pair_level <- any(n_levels > 1L)
        if( pair_level ){
            message("'split.by' varies within subjects: each pair takes ",
                "the value of its earlier sample.")
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
            mid = (time[k] + time[l]) / 2, d = d[ut], k = k, l = l)
    })
    pairs <- do.call(rbind, unname(res))
    pairs$d <- .dissimilarity_scale(pairs$d, squared, bounded)
    if( !is.null(split.by) ){
        # The level of the earlier sample
        first <- ifelse(time[pairs$k] <= time[pairs$l], pairs$k, pairs$l)
        pairs$level <- level[first]
        attr(pairs, "pair.level") <- pair_level
    }
    # Columns of the input that k and l index
    attr(pairs, "samples") <- samples
    return(pairs)
}

# The pairs whose earlier sample is TRUE in the logical column subset.by
.subset_pairs <- function(pairs, x, time.col, subset.by){
    samples <- attr(pairs, "samples")
    pair_level <- attr(pairs, "pair.level")
    keep <- x[[subset.by]][samples]
    time <- x[[time.col]][samples]
    first <- ifelse(time[pairs$k] <= time[pairs$l], pairs$k, pairs$l)
    pairs <- pairs[keep[first] %in% TRUE, , drop = FALSE]
    if( nrow(pairs) == 0L ){
        stop("No pair starts from a sample with 'subset.by' TRUE.",
            call. = FALSE)
    }
    attr(pairs, "samples") <- samples
    attr(pairs, "pair.level") <- pair_level
    return(pairs)
}

# Squared, or on the scale -log(1 - d) for bounded dissimilarities
.dissimilarity_scale <- function(d, squared, bounded){
    if( squared ){
        d <- d^2
    }
    if( bounded ){
        if( any(d < 0 | d > 1, na.rm = TRUE) ){
            stop("'bounded = TRUE' needs dissimilarities between 0 and 1.",
                call. = FALSE)
        }
        d <- -log1p(-pmin(d, 1 - 1e-9))
    }
    return(d)
}

# Mean dissimilarity between samples of different subjects, every pair of
# subjects counting equally
.between_subject_level <- function(
        x, time.col, group, assay.type, method, split.by, squared, bounded,
        args){
    keep <- !is.na(x[[time.col]]) & !is.na(x[[group]])
    if( !is.null(split.by) ){
        keep <- keep & !is.na(x[[split.by]])
    }
    x <- x[, keep]
    subject <- as.character(x[[group]])
    if( length(unique(subject)) < 2L ){
        stop("'relative = TRUE' needs at least two subjects.", call. = FALSE)
    }
    d <- do.call(getDissimilarity, c(
        list(x = t(assay(x, assay.type)), method = method), args))
    d <- .dissimilarity_scale(as.matrix(d), squared, bounded)
    # Sums and numbers of sample pairs for every pair of subjects
    sums <- rowsum(t(rowsum(d, subject)), subject)
    n <- as.numeric(table(subject)[rownames(sums)])
    means <- sums / outer(n, n)
    B <- mean(means[upper.tri(means)])
    return(B)
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
        pairs, lags, split.by, time.point, window, conf.level,
        permutations = 0){
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
    if( isTRUE(attr(pairs, "pair.level")) ){
        comp <- lapply(seq_along(lev)[-1L], function(i){
            .compare_pair_levels(
                pairs, lev[[1L]], lev[[i]], lags, window, a0, est[[1L]],
                est[[i]])
        })
    } else{
        comp <- lapply(est[-1L], function(e2){
            .compare_levels(est[[1L]], e2)
        })
    }
    if( permutations > 0 ){
        for( i in seq_along(comp) ){
            comp[[i]]$p.value.permutation <- .permutation_p(
                pairs, lev[[1L]], lev[[i + 1L]], lags, window, conf.level,
                a0, comp[[i]], permutations)
        }
        est <- lapply(est, function(e){
            e$p.value.permutation <- NA_real_
            e
        })
    }
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
    f <- .local_design(pairs, l, a0, k, window)
    res <- data.frame(
        lag = l, estimate = NA_real_, se = NA_real_, df = NA_real_,
        n.subjects = f$G, n.pairs = length(f$rows), supported = FALSE,
        time.point = a0)
    if( !f$ok ){
        return(res)
    }
    G <- f$G
    d <- pairs$d[f$rows]
    b <- f$Ai %*% crossprod(f$X, f$w * d)
    r <- d - drop(f$X %*% b)
    U <- rowsum(f$X * (f$w * r), f$subj)
    V <- f$Ai %*% crossprod(U) %*% f$Ai * G / (G - 1)
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

# The pairs in the window of lag l (and time a0), the support rules, and the
# design of the local regression
.local_design <- function(pairs, l, a0, k, window){
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
    res <- list(ok = FALSE, rows = rows, G = G)
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
    # The subjects' scores sum to zero, so the sandwich needs clearly more
    # subjects than coefficients
    if( G < 2L * ncol(X) ){
        return(res)
    }
    # Each subject's pairs in the window count as one subject
    w <- as.numeric(1 / ns[subj])
    A <- crossprod(X, X * w)
    if( qr(A)$rank < ncol(X) ){
        return(res)
    }
    res <- list(
        ok = TRUE, rows = rows, G = G, X = X, w = w, subj = subj,
        Ai = solve(A))
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

# Ratio of D(L) between two levels of a sample-level state: one regression on
# the pairs of both levels, with an indicator of the second level, lag terms
# for each level and time terms shared by the levels. The log ratio is
# log((b0 + b1) / b0), with a cluster-robust (CR1) standard error by the
# delta method.
.compare_pair_levels <- function(pairs, lev1, lev2, lags, window, a0, e1, e2){
    p <- pairs[pairs$level %in% c(lev1, lev2), , drop = FALSE]
    res <- lapply(seq_along(lags), function(i){
        l <- lags[[i]]
        out <- data.frame(
            lag = l, estimate = NA_real_, se = NA_real_, df = NA_real_,
            n.subjects = NA_integer_, n.pairs = NA_integer_,
            supported = FALSE, time.point = a0[[i]], p.value = NA_real_,
            type = "ratio")
        f <- .pair_level_design(p, lev2, l, a0[[i]], window)
        G <- f$G
        out$n.subjects <- G
        out$n.pairs <- length(f$rows)
        # Each level must be supported on its own
        if( !e1$supported[[i]] || !e2$supported[[i]] || !f$ok ){
            return(out)
        }
        d <- p$d[f$rows]
        b <- drop(f$Ai %*% crossprod(f$X, f$w * d))
        if( b[[1L]] <= 0 || b[[1L]] + b[[2L]] <= 0 ){
            return(out)
        }
        r <- d - drop(f$X %*% b)
        U <- rowsum(f$X * (f$w * r), f$subj)
        V <- f$Ai %*% crossprod(U) %*% f$Ai * G / (G - 1)
        lr <- log(b[[1L]] + b[[2L]]) - log(b[[1L]])
        grad <- c(1 / (b[[1L]] + b[[2L]]) - 1 / b[[1L]],
            1 / (b[[1L]] + b[[2L]]), rep(0, ncol(f$X) - 2L))
        se <- sqrt(drop(t(grad) %*% V %*% grad))
        out$estimate <- exp(lr)
        out$se <- se
        out$df <- G - 1
        out$supported <- TRUE
        out$p.value <- 2 * stats::pt(-abs(lr / se), G - 1)
        return(out)
    })
    res <- do.call(rbind, res)
    return(res)
}

# Design of the regression on the pairs of two levels of a sample-level state
# in the window of lag l (and time a0); the second coefficient is that of the
# indicator of lev2
.pair_level_design <- function(p, lev2, l, a0, window){
    rows <- which(p$h >= l / window & p$h <= l * window)
    if( !is.na(a0) ){
        r <- (1 + pmax(p$mid[rows], 0)) / (1 + max(a0, 0))
        rows <- rows[r >= 1 / 2 & r <= 2]
    }
    q <- p[rows, , drop = FALSE]
    subj <- q$subject
    ns <- table(subj)
    G <- length(ns)
    res <- list(ok = FALSE, rows = rows, G = G)
    g <- as.numeric(q$level == lev2)
    u <- log(q$h / l)
    X <- cbind(1, g, u, u^2, g * u, g * u^2)
    if( !is.na(a0) ){
        z <- log1p(pmax(q$mid, 0)) - log1p(max(a0, 0))
        X <- cbind(X, z, z^2)
    }
    if( G < 2L * ncol(X) ){
        return(res)
    }
    w <- as.numeric(1 / ns[subj])
    A <- crossprod(X, X * w)
    if( qr(A)$rank < ncol(X) ){
        return(res)
    }
    res <- list(
        ok = TRUE, rows = rows, G = G, X = X, w = w, subj = subj,
        Ai = solve(A))
    return(res)
}

# Permutation p-value of each lag's comparison between two levels: the
# levels are reassigned among the subjects of the two levels, keeping their
# numbers, and the t statistic of the log ratio is recomputed
.permutation_p <- function(
        pairs, lev1, lev2, lags, window, conf.level, a0, obs,
        permutations){
    p <- pairs[pairs$level %in% c(lev1, lev2), , drop = FALSE]
    subj <- unique(data.frame(
        subject = p$subject, level = as.character(p$level)))
    t_obs <- log(obs$estimate) / obs$se
    t_perm <- vapply(seq_len(permutations), function(b){
        lab <- stats::setNames(sample(subj$level), subj$subject)
        pl <- lab[p$subject]
        e1 <- .local_variogram(
            p[pl == lev1, , drop = FALSE], lags, window, conf.level, a0)
        e2 <- .local_variogram(
            p[pl == lev2, , drop = FALSE], lags, window, conf.level, a0)
        cmp <- .compare_levels(e1, e2)
        ifelse(cmp$supported, log(cmp$estimate) / cmp$se, NA_real_)
    }, numeric(length(lags)))
    t_perm <- matrix(t_perm, nrow = length(lags))
    res <- vapply(seq_along(lags), function(i){
        v <- t_perm[i, ]
        v <- v[!is.na(v)]
        if( !obs$supported[[i]] || length(v) == 0L ){
            return(NA_real_)
        }
        (1 + sum(abs(v) >= abs(t_obs[[i]]))) / (1 + length(v))
    }, numeric(1L))
    return(res)
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
        "p.value", intersect("p.value.permutation", names(res)),
        "n.subjects", "n.pairs", "supported")
    if( !is.null(time.point) ){
        cols <- c(cols, "time.point")
    }
    res <- res[, cols]
    return(res)
}
