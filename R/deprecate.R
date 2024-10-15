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

#' @rdname deprecate
#' @export
setGeneric("getStepwiseDivergence", signature = c("x"), function(x, ... )
    standardGeneric("getStepwiseDivergence"))

#' @rdname deprecate
#' @export
setMethod("getStepwiseDivergence", signature = c(x = "ANY"), function(x, ...){
    .Deprecated(msg = "'getStepwiseDivergence' is deprecated. 
                Use 'addStepwiseDivergence' instead.")
    addStepwiseDivergence(x, ...)
}
)

#' @rdname deprecate
#' @export
setGeneric("getBaselineDivergence", signature = c("x"), function(x, ... )
    standardGeneric("getBaselineDivergence"))

#' @rdname deprecate
#' @export
setMethod("getBaselineDivergence", signature = c(x = "ANY"), function(x, ...){
    .Deprecated(msg = "'getBaselineDivergence' is deprecated. 
                Use 'addBaselineDivergence' instead.")
    addBaselineDivergence(x, ...)
}
)