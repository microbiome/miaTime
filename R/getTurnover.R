#' @name
#' getTurnover
#'
#' @export
#'
#' @title
#' Turnover of features between consecutive samples
#'
#' @description
#' Calculates the proportion of features that appear or disappear between each
#' sample and the previous sample of the same subject.
#'
#' @details
#' Turnover (Hallett et al. 2016) compares which features are present (above
#' \code{detection}) in two samples. With \eqn{G} features gained, \eqn{L}
#' features lost and \eqn{N} features present in either sample, total turnover
#' is \eqn{(G + L) / N}, appearance is \eqn{G / N} and disappearance is
#' \eqn{L / N}. The previous sample is chosen as in
#' \code{\link[miaTime:getStepwiseDivergence]{getStepwiseDivergence()}}.
#'
#' @references
#' Hallett LM, et al. (2016) codyn: An R package of community dynamics
#' metrics. Methods in Ecology and Evolution 7(10), 1146-1151.
#' doi: 10.1111/2041-210X.12569
#'
#' @return
#' \code{getTurnover} returns a \code{DataFrame} with one row per sample and
#' the turnover, the time from the previous sample and the previous sample.
#' \code{addTurnover} adds these to \code{colData(x)}. The first sample of each
#' subject gets \code{NA}.
#'
#' @inheritParams getStepwiseDivergence
#'
#' @param method \code{Character scalar}. The turnover to calculate:
#' \code{"total"}, \code{"appearance"} or \code{"disappearance"}.
#' (Default: \code{"total"})
#'
#' @param detection \code{Numeric scalar}. Abundance above which a feature is
#' considered present. (Default: \code{0})
#'
#' @param name \code{Character vector}. Specifies the names of the three
#' columns of results. (Default: \code{c("turnover", "time_diff",
#' "ref_samples")})
#'
#' @param ... additional arguments.
#'
#' @examples
#' library(miaTime)
#'
#' # Two faecal samples of each subject in an IBD cohort
#' data(HMP_2019_ibdmdb)
#' tse <- HMP_2019_ibdmdb
#'
#' # Proportion of taxa gained since the first sample
#' tse <- addTurnover(tse, time.col = "days_from_first_collection",
#'     group = "subject_id", assay.type = "relabundance",
#'     method = "appearance", name = c("gained", "time_diff", "ref_samples"))
#'
#' # Proportion of taxa lost
#' tse[["lost"]] <- getTurnover(tse, time.col = "days_from_first_collection",
#'     group = "subject_id", assay.type = "relabundance",
#'     method = "disappearance")[["turnover"]]
#'
#' df <- as.data.frame(colData(tse))
#' df <- df[!is.na(df[["gained"]]), ]
#' head(df[, c("subject_id", "disease_subtype", "time_diff", "gained", "lost")])
#'
#' @seealso
#' \code{\link[miaTime:getStepwiseDivergence]{getStepwiseDivergence()}}
NULL

#' @rdname getTurnover
#' @export
setMethod("getTurnover", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, assay.type = "counts", time.interval = 1L,
            time.lag = NULL, group = NULL,
            method = c("total", "appearance", "disappearance"),
            detection = 0, ...){
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
        method <- match.arg(method)
        temp <- .check_input(detection, list("numeric scalar"))
        if( is.null(colnames(x)) ){
            colnames(x) <- paste0("col", seq_len(ncol(x)))
        }
        ########################### INPUT CHECK END ############################
        x <- .check_and_get_altExp(x, ...)
        args <- .add_reference_samples_to_coldata(
            x, time.col, group, time.interval = time.interval,
            time.lag = time.lag, reference.method = "stepwise", ...)
        present <- assay(args[["x"]], assay.type) > detection
        refs <- args[["x"]][[args[["reference"]]]]
        # With repeated time points, a sample has several previous samples;
        # their turnovers are averaged
        res <- vapply(seq_len(ncol(present)), function(j){
            r <- unlist(refs[j])
            r <- r[!is.na(r)]
            if( length(r) == 0L ){
                return(NA_real_)
            }
            mean(vapply(r, function(k)
                .turnover(present[, j], present[, k], method), numeric(1)))
        }, numeric(1))
        time_res <- .get_time_difference(
            args[["x"]], time.col, args[["reference"]])
        res <- .convert_divergence_to_df(
            args[["x"]], res, time_res, args[["reference"]],
            name = c("turnover", "time_diff", "ref_samples"))
        return(res)
    }
)

#' @rdname getTurnover
#' @export
setMethod("addTurnover", signature = c(x = "SummarizedExperiment"),
    function(x, name = c("turnover", "time_diff", "ref_samples"), ...){
        .check_input(name, c("character vector"), length = 3L)
        res <- getTurnover(x, ...)
        res <- as.list(res) |> unname()
        x <- .add_values_to_colData(x, res, name, ...)
        return(x)
    }
)

################################ HELP FUNCTIONS ################################

# Proportion of the features present in either sample that were gained, lost
# or both, from the previous sample (prev) to the current one (cur).
.turnover <- function(cur, prev, method){
    n <- sum(cur | prev)
    if( n == 0L ){
        return(NA_real_)
    }
    changed <- switch(method,
        total = cur != prev,
        appearance = cur & !prev,
        disappearance = !cur & prev)
    return(sum(changed) / n)
}
