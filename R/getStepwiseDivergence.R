#' Beta diversity between consecutive time steps 
#'
#' Calculates sample dissimilarity between consecutive time points (t, t+i),
#' within a group (subject, reaction chamber, or similar). The corresponding
#' time difference is returned as well. The method operates on
#' `SummarizedExperiment` objects, and the results are stored in `colData`.
#'
#' @inheritParams addBaselineDivergence
#' 
#' @param time.interval \code{Integer scalar}. Indicates the increment between 
#' time steps. If you need to take every second, every third, or so, time step 
#' only, then increase this accordingly. (Default: \code{1L})
#' 
#' @return a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' containing the sample dissimilarity and corresponding time difference between
#' samples (across n time steps), within each level of the grouping factor.
#'
#' @name addStepwiseDivergence
#' @export
#'
#' @examples
#' library(miaTime)
#'
#' data(hitchip1006)
#' tse <- transformAssay(hitchip1006, method = "relabundance")
#' 
#' # Calculate divergence
#' tse <- addStepwiseDivergence(
#'     tse, group = "subject",
#'     time_interval = 1,
#'     time.col = "time",
#'     assay.type="relabundance"
#'     )
#' 
NULL

#' @rdname addStepwiseDivergence
#' @export
#'
setGeneric("getStepwiseDivergence", signature = c("x"), function(x, ...)
  standardGeneric("getStepwiseDivergence"))

#' @rdname addStepwiseDivergence
#' @export
setMethod("getStepwiseDivergence", signature = c(x = "ANY"),
    function(
        x,
        time.col,
        assay.type = "counts",
        time.interval = 1L,
        group = NULL,
        method = "bray",
        name = "divergence",
        name.time = "time_diff",
        ...){
        ############################# INPUT CHECK ##############################
        temp <- .check_input(
            time.col, list("character scalar"), colnames(colData(x)))
        if( !is.numeric(x[[time.col]]) ){
            stop("'time.col' must specify numeric column from colData(x)",
                call. = FALSE)
        }
        #
        .check_assay_present(assay.type, x)
        #
        temp <- .check_input(
            group, list(NULL, "character scalar"), colnames(colData(x)))
        #
        temp <- .check_input(method, list("character scalar"))
        #
        temp <- .check_input(name, list(NULL, "character scalar"))
        #
        temp <- .check_input(name.time, list(NULL, "character scalar"))
        #
        if( is.null(rownames(x)) ){
            rownames(x) <- paste0("row", seq_len(nrow(x)))
        }
        if( is.null(colnames(x)) ){
            colnames(x) <- paste0("col", seq_len(ncol(x)))
        }
        ########################### INPUT CHECK END ############################
        x <- .check_and_get_altExp(x, ...)
        # Add stepwise samples to colData
        x <- .add_reference_samples_to_coldata(
            x, time.col, group, time.interval = time.interval,
            reference.method = "stepwise", ...)
        reference <- x[[2]]
        x <- x[[1]]
        # Calculate divergences
        res <- getDivergence(
            x, assay.type = assay.type, reference = reference,
            method = method, ...)
        # Add time difference
        time_res <- .get_time_difference(x, time.col, reference)
        # Create a DF to return
        res <- .convert_divergence_to_df(x, res, time_res, name, name.time)
        return(res)
    }
)

#' @rdname addStepwiseDivergence
#' @export
setGeneric("addStepwiseDivergence", signature = "x", function(x, ...)
    standardGeneric("addStepwiseDivergence"))

#' @rdname addStepwiseDivergence
#' @export
setMethod("addStepwiseDivergence", signature = c(x = "SummarizedExperiment"),
    function(x, name = "divergence", name.time = "time_diff", ...){
      # Calculate divergence
        res <- getStepwiseDivergence(x,  ...)
        # Add to colData
        res <- as.list(res) |> unname()
        x <- .add_values_to_colData(x, res, list(name, name.time), ...)
        return(x)
    }
)
