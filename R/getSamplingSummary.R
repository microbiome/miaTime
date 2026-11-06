#' @name
#' getSamplingSummary
#'
#' @title
#' Summarize sampling over time
#'
#' @description
#' Summarizes how each group, for instance a subject, was sampled over time.
#'
#' @details
#' Temporal metrics, such as stepwise divergence or stability, depend on how
#' many time points each group has and how far apart they are. This function
#' reports, for each group, the number of samples and distinct time points,
#' the first and last time point, the time span, and the minimum, median and
#' maximum gap between consecutive time points. Gaps are \code{NA} for groups
#' with a single time point.
#'
#' @return
#' \code{DataFrame} with one row per group.
#'
#' @inheritParams addBaselineDivergence
#'
#' @param group \code{NULL} or \code{Character scalar}. Specifies a column
#' from \code{colData(x)} that defines the groups. If \code{NULL}, all
#' samples are summarized together. (Default: \code{NULL})
#'
#' @param ... additional arguments. Currently not used.
#'
#' @examples
#' library(miaTime)
#'
#' data(hitchip1006)
#' res <- getSamplingSummary(hitchip1006, time.col = "time", group = "subject")
#'
#' # Number of subjects by the number of time points
#' table(res$n_time_points)
#'
NULL

#' @rdname getSamplingSummary
#' @export
setMethod("getSamplingSummary", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group = NULL, ...){
        ############################# INPUT CHECK ##############################
        temp <- .check_input(
            time.col, list("character scalar"), colnames(colData(x)))
        if( !is.numeric(x[[time.col]]) ){
            stop("'time.col' must specify numeric column from colData(x)",
                call. = FALSE)
        }
        temp <- .check_input(
            group, list(NULL, "character scalar"), colnames(colData(x)))
        ########################### INPUT CHECK END ############################
        time <- x[[time.col]]
        groups <- if( is.null(group) ) rep("all", ncol(x)) else x[[group]]
        res <- lapply(split(time, groups), .summarize_sampling)
        res <- do.call(rbind, res) |> DataFrame()
        return(res)
    }
)

################################ HELP FUNCTIONS ################################

# Summarize the time points of a single group
#' @importFrom stats median
.summarize_sampling <- function(time){
    time_points <- sort(unique(time[!is.na(time)]))
    gaps <- diff(time_points)
    res <- data.frame(
        n_samples = length(time),
        n_time_points = length(time_points),
        start = if( length(time_points) ) min(time_points) else NA_real_,
        end = if( length(time_points) ) max(time_points) else NA_real_,
        min_gap = if( length(gaps) ) min(gaps) else NA_real_,
        median_gap = if( length(gaps) ) median(gaps) else NA_real_,
        max_gap = if( length(gaps) ) max(gaps) else NA_real_
    )
    res[["span"]] <- res[["end"]] - res[["start"]]
    res <- res[, c("n_samples", "n_time_points", "start", "end", "span",
        "min_gap", "median_gap", "max_gap")]
    return(res)
}
