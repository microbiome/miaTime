#' @title Get Time Series Measures
#' @description Calculate various time series measures (e.g., ACF, PACF) 
#' from a specified column in the `colData` of a `TreeSE` object.
#' 
#' @param x A
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#' @param time.field \code{Character scalar}. Specifies the name of the
#' time series field in `colData`. (Default: \code{"Time.hr"})
#' @param measures \code{Character vector}. Specifies the time series measures 
#' to compute. Supported measures include "acf", "pacf","Box.test", and "arima".
#' (Default: \code{ "c("acf", "pacf")"})
#' @param ... Additional arguments passed to specific time series functions.
#' 
#' @return A list containing the results of the requested time series measures.
#' 
#' @examples
#' 
#' data(minimalgut)
#' tse <- minimalgut
#' result <- getTsMeasure(tse, time.field = "Time.hr", measure = "acf")
#' 
#' @rdname getTsMeasure
#' @export
setGeneric("getTsMeasure", signature = "x",
    function(x, ...) 
       standardGeneric("getTsMeasure")
)

#' @rdname getTsMeasure
#' @export
setMethod("getTsMeasure", signature = c(x = "TreeSummarizedExperiment"),
    function(x, time.field = "Time.hr", measures = c("acf", "pacf"), ...) {
        ############################## Input check #############################
        # Check if the specified time.field exists in colData
        if(!is.character(time.field)){
            stop("'time.field' must be a on empty single character value.", call. = FALSE)
        }
        if (!(time.field %in% colnames(colData(x)))) {
            stop("'time.field' does not exist in colData of the TreeSE 
                object.", call. = FALSE)
        }
        
        ############################ Input check end ###########################
        
        ############################ Measure Calculation ######################
        results <- lapply(measures, function(measure) {
            .getTsMeasures(x, time.field, measure, ...)
        })
        names(results) <- measures
        ############################ Measure Calculation end ###################
        
        return(results)
        }
)
