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
#' (Default: \code{ ""acf""})
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
    function(x, time.field = "Time.hr", measures = "acf", ...) {
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
            .get_ts_measures(x, time.field, measure, ...)
        })
        names(results) <- measures
        ############################ Measure Calculation end ###################
        
        return(results)
        }
)

#' @importFrom stats acf pacf Box.test arima
# wrapper for time series measures
.get_ts_measures <- function(x, time.field, measure, ...) {
    # Extract the time series data from the specified column in colData
    ts_data <- colData(x)[[time.field]]
    
    # Ensure time series is numeric
    if (!is.numeric(ts_data)) {
        stop("'time.field' does not contain numeric time series 
                data.", call. = FALSE)
    }
    
    # Get the correct ts measure to call
    FUN <- switch(measure,
                  "acf" = function() acf(ts_data, plot = FALSE),
                  "pacf" = function() pacf(ts_data, plot = FALSE),
                  "Box.test" = function() Box.test(ts_data, type = "Ljung-Box"),
                  "arima" = function() arima(ts_data),
                  stop("Unsupported measure: ", measure, call. = FALSE))
    
    # Calculate the result
    res <- FUN()
    
    # Process the result based on the measure
    if (measure %in% c("acf", "pacf")) {
        res <- res$acf
    }
    
    return(res)
}
