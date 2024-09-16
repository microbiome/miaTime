#' @rdname deprecate
#' @export
setGeneric("getTimeDivergence", signature = c("x"), function(x, ... )
    standardGeneric("getTimeDivergence"))

#' @rdname deprecate
#' @export
setMethod("getTimeDivergence", signature = c(x = "ANY"), function(x, ...){
    .Deprecated(msg = "test here")
    addStepwiseDivergence(x, ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric("getStepwiseDivergence", signature = c("x"), function(x, ... )
    standardGeneric("getTimeDivergence"))

#' @rdname deprecate
#' @export
setMethod("getStepwiseDivergence", signature = c(x = "ANY"), function(x, ...){
    .Deprecated(msg = "Text here")
    addStepwiseDivergence(x, ...)
}
)

#' @rdname deprecate
#' @export
setGeneric("getBaselineDivergence", signature = c("x"), function(x, ... )
    standardGeneric("getTimeDivergence"))

#' @rdname deprecate
#' @export
setMethod("getTimeDivergence", signature = c(x = "ANY"), function(x, ...){
    .Deprecated(msg = "add rexxr")
    addBaselineDivergence(x, ...)
}
)