#' @name
#' getRankShift
#'
#' @export
#'
#' @title
#' Mean rank shift between consecutive samples
#'
#' @description
#' Calculates how much the abundance ranks of features change between each
#' sample and the previous sample of the same subject.
#'
#' @details
#' The mean rank shift (Collins et al. 2008) compares two samples through the
#' features present (above zero) in both. The features are ranked by
#' abundance in each sample, with tied values getting their average rank, and
#' the rank shift is the mean absolute change in rank. It is 0 when the order
#' of the features does not change. The previous sample is chosen as in
#' \code{\link[miaTime:getStepwiseDivergence]{getStepwiseDivergence()}}.
#'
#' @references
#' Collins SL, et al. (2008) Rank clocks and plant community dynamics.
#' Ecology 89(12), 3534-3541. doi: 10.1890/07-1646.1
#'
#' @return
#' \code{getRankShift} returns a \code{DataFrame} with one row per sample and
#' the rank shift, the time from the previous sample and the previous
#' sample. \code{addRankShift} adds these to \code{colData(x)}. The first
#' sample of each subject gets \code{NA}.
#'
#' @inheritParams getStepwiseDivergence
#'
#' @param name \code{Character vector}. Specifies the names of the three
#' columns of results. (Default: \code{c("rank_shift", "time_diff",
#' "ref_samples")})
#'
#' @param ... additional arguments.
#'
#' @examples
#' library(miaTime)
#'
#' # A defined community of 16 species in three replicate bioreactors
#' data(minimalgut)
#' tse <- addRankShift(minimalgut, time.col = "Time_hr_num",
#'     group = "StudyIdentifier")
#'
#' # Rank shift over time in each bioreactor; the first samples have no
#' # previous sample
#' library(ggplot2)
#' df <- as.data.frame(colData(tse))
#' df <- df[!is.na(df[["rank_shift"]]), ]
#' ggplot(df, aes(x = Time_hr_num, y = rank_shift, colour = StudyIdentifier)) +
#'     geom_line()
#'
#' @seealso
#' \code{\link[miaTime:getStepwiseDivergence]{getStepwiseDivergence()}}
NULL

#' @rdname getRankShift
#' @export
setMethod("getRankShift", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, assay.type = "counts", time.interval = 1L,
            time.lag = NULL, group = NULL, ...){
        ############################# INPUT CHECK ##############################
        temp <- .check_input(
            time.col, list("character scalar"), colnames(colData(x)))
        if( !is.numeric(x[[time.col]]) ){
            stop("'time.col' must specify numeric column from colData(x)",
                call. = FALSE)
        }
        .check_assay_present(assay.type, x)
        temp <- .check_input(
            group, list(NULL, "character scalar"), colnames(colData(x)))
        if( is.null(colnames(x)) ){
            colnames(x) <- paste0("col", seq_len(ncol(x)))
        }
        ########################### INPUT CHECK END ############################
        x <- .check_and_get_altExp(x, ...)
        args <- .add_reference_samples_to_coldata(
            x, time.col, group, time.interval = time.interval,
            time.lag = time.lag, reference.method = "stepwise", ...)
        mat <- assay(args[["x"]], assay.type)
        refs <- args[["x"]][[args[["reference"]]]]
        # With repeated time points, a sample has several previous samples;
        # their rank shifts are averaged
        res <- vapply(seq_len(ncol(mat)), function(j){
            r <- unlist(refs[j])
            r <- r[!is.na(r)]
            if( length(r) == 0L ){
                return(NA_real_)
            }
            mean(vapply(r, function(k) .rank_shift(mat[, j], mat[, k]),
                numeric(1)))
        }, numeric(1))
        time_res <- .get_time_difference(
            args[["x"]], time.col, args[["reference"]])
        res <- .convert_divergence_to_df(
            args[["x"]], res, time_res, args[["reference"]],
            name = c("rank_shift", "time_diff", "ref_samples"))
        return(res)
    }
)

#' @rdname getRankShift
#' @export
setMethod("addRankShift", signature = c(x = "SummarizedExperiment"),
    function(x, name = c("rank_shift", "time_diff", "ref_samples"), ...){
        .check_input(name, c("character vector"), length = 3L)
        res <- getRankShift(x, ...)
        res <- as.list(res) |> unname()
        x <- .add_values_to_colData(x, res, name, ...)
        return(x)
    }
)

################################ HELP FUNCTIONS ################################

# Mean absolute change in the abundance rank of the features present in both
# samples.
.rank_shift <- function(a, b){
    present <- a > 0 & b > 0
    if( !any(present) ){
        return(NA_real_)
    }
    return(mean(abs(rank(a[present]) - rank(b[present]))))
}
