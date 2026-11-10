#' @name
#' getAutocorrelation
#'
#' @export
#'
#' @title
#' Autocorrelation of features over time
#'
#' @description
#' Calculates the autocorrelation of each feature within each subject at
#' lags of one or more time steps.
#'
#' @details
#' The samples of each subject must lie on a regular grid of time steps. The
#' step is the shortest time between two samples of a subject. Grid points
#' without a sample are treated as missing, and replicated samples at the
#' same time are averaged. The autocorrelation is computed with
#' \code{\link[stats:acf]{acf}}, passing missing values through.
#'
#' For irregularly sampled data, where samples do not lie on a common grid,
#' use \code{getVariogram()}, which estimates how
#' fast communities change at arbitrary time lags.
#'
#' @return
#' \code{getAutocorrelation} returns a \code{DataFrame} with one row per
#' feature, subject and lag, and the columns \code{FeatureID}, \code{group}
#' (when given), \code{lag} (in time steps), \code{time.lag} (in the units of
#' \code{time.col}), \code{acf} and \code{n.pairs} (the pairs of samples that
#' contribute to the estimate). \code{addAutocorrelation} returns \code{x}
#' with the table in its \code{metadata}.
#'
#' @inheritParams addBaselineDivergence
#'
#' @param max.lag \code{NULL} or \code{Integer scalar}. The largest lag, in
#' time steps. If \code{NULL}, the default of \code{acf}.
#' (Default: \code{NULL})
#'
#' @param name \code{Character scalar}. Specifies a name for storing the
#' results. (Default: \code{"autocorrelation"})
#'
#' @param ... additional arguments.
#'
#' @examples
#' library(miaTime)
#'
#' # Artificial gut with four vessels sampled every hour
#' data(SilvermanAGutData)
#' tse <- SilvermanAGutData
#' tse <- tse[, tse$SampleType == "Hourly" &
#'     tse$Normal_Noise_Sample == "Normal"]
#' # Samples named by day and hour, with reads
#' tse <- tse[, grepl("^T\\.+[0-9]+d\\.[0-9]+h", colnames(tse)) &
#'     colSums(assay(tse, "counts")) > 0]
#' # Hours since inoculation, from the sample names
#' days <- as.numeric(sub("^T\\.+([0-9]+)d.*", "\\1", colnames(tse)))
#' hours <- as.numeric(sub("^T\\.+[0-9]+d\\.([0-9]+)h.*", "\\1", colnames(tse)))
#' tse$hour <- days * 24 + hours
#'
#' tse <- transformAssay(tse, method = "relabundance")
#' df <- getAutocorrelation(tse, time.col = "hour", group = "Vessel",
#'     assay.type = "relabundance", max.lag = 24)
#'
#' # Mean autocorrelation of the five most abundant taxa over the vessels
#' top <- getTop(tse, top = 5, assay.type = "relabundance")
#' df <- df[df[["FeatureID"]] %in% top, ]
#' aggregate(acf ~ FeatureID + time.lag, data = as.data.frame(df),
#'     FUN = mean)[1:10, ]
#'
#' @seealso
#' \code{\link[stats:acf]{acf}}
NULL

#' @rdname getAutocorrelation
#' @export
setMethod("addAutocorrelation", signature = c(x = "SummarizedExperiment"),
    function(x, name = "autocorrelation", ...){
        temp <- .check_input(name, list("character scalar"))
        x <- .check_and_get_altExp(x, ...)
        args <- c(list(x = x), list(...)[!names(list(...)) %in% c("altexp")])
        res <- do.call(getAutocorrelation, args)
        x <- .add_values_to_metadata(x, name, res, ...)
        return(x)
    }
)

#' @rdname getAutocorrelation
#' @export
setMethod("getAutocorrelation", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, assay.type = "counts", group = NULL,
            max.lag = NULL, ...){
        ############################## Input check #############################
        x <- .check_and_get_altExp(x, ...)
        temp <- .check_input(
            time.col, list("character scalar"), colnames(colData(x)))
        if( !is.numeric(x[[time.col]]) || anyNA(x[[time.col]]) ){
            stop("'time.col' must specify a numeric column from colData(x) ",
                "without missing values.", call. = FALSE)
        }
        .check_assay_present(assay.type, x)
        temp <- .check_input(
            group, list(NULL, "character scalar"), colnames(colData(x)))
        temp <- .check_input(max.lag, list(NULL, "numeric scalar"))
        ########################### Input check end ############################
        time <- x[[time.col]]
        subject <- if( is.null(group) ) rep("all", ncol(x)) else
            as.character(x[[group]])
        mat <- as.matrix(assay(x, assay.type))
        res <- .calculate_autocorrelation(mat, time, subject, max.lag)
        if( is.null(group) ){
            res[["group"]] <- NULL
        } else{
            names(res)[names(res) == "group"] <- group
        }
        return(DataFrame(res, check.names = FALSE))
    }
)

################################ HELP FUNCTIONS ################################

# Autocorrelation of each feature (row of mat) within each subject.
#' @importFrom stats acf na.pass
.calculate_autocorrelation <- function(mat, time, subject, max.lag){
    # The step is the shortest time between two samples of a subject.
    gaps <- unlist(tapply(time, subject, function(t) diff(sort(unique(t)))))
    step <- min(gaps[gaps > 0])
    # Replicated samples at the same time are averaged.
    key <- paste(subject, time, sep = "\r")
    if( anyDuplicated(key) ){
        message("The dataset contains replicated samples. The average ",
            "abundance is calculated for each time point.")
        mat <- t(rowsum(t(mat), key, reorder = FALSE)) /
            rep(as.vector(table(factor(key, levels = unique(key)))),
                each = nrow(mat))
        first <- !duplicated(key)
        time <- time[first]
        subject <- subject[first]
    }
    res <- lapply(split(seq_along(time), subject), function(idx){
        k <- (time[idx] - min(time[idx])) / step
        if( any(abs(k - round(k)) > 1e-8) ){
            stop("The samples of each subject must lie on a regular grid of ",
                "time steps; use getVariogram() for irregular sampling.",
                call. = FALSE)
        }
        k <- round(k)
        series <- matrix(NA_real_, nrow(mat), max(k) + 1L)
        series[, k + 1L] <- mat[, idx, drop = FALSE]
        n <- ncol(series)
        lags <- seq_len(min(max.lag %||% floor(10 * log10(n)), n - 1L))
        acfs <- vapply(seq_len(nrow(series)), function(i){
            acf(series[i, ], lag.max = max(lags), plot = FALSE,
                na.action = na.pass)[["acf"]][lags + 1L]
        }, numeric(length(lags)))
        present <- !is.na(series[1L, ])
        n_pairs <- vapply(lags, function(l){
            sum(present[seq_len(n - l)] & present[(l + 1L):n])
        }, numeric(1))
        data.frame(
            FeatureID = rep(rownames(mat), each = length(lags)),
            group = subject[idx[1L]],
            lag = rep(lags, nrow(mat)),
            time.lag = rep(lags * step, nrow(mat)),
            acf = as.vector(acfs),
            n.pairs = rep(n_pairs, nrow(mat)))
    })
    do.call(rbind, unname(res))
}
