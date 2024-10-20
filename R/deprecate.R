#' These functions are deprecated. Please use other functions instead.
#' 
#' @param x A
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#' 
#' @param ... Additional parameters. See dedicated function.
#' 
#' @name deprecate
NULL

#' @rdname deprecate
#' @export
setGeneric("getTimeDivergence", signature = c("x"), function(x, ... )
    standardGeneric("getTimeDivergence"))

#' @rdname deprecate
#' @export
setMethod("getTimeDivergence", signature = c(x = "ANY"), function(x, ...){
    .Deprecated(msg = "'getTimeDivergence' is deprecated. 
        Use 'addStepwiseDivergence' instead.")
    addStepwiseDivergence(x, ...)
    }
)

