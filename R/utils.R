################################################################################
# internal methods loaded from other packages

.check_altExp_present <- mia:::.check_altExp_present
.calc_reference_dist <- mia:::.calc_reference_dist
.get_mat_from_sce <- scater:::.get_mat_from_sce

################################################################################
# internal wrappers for getter/setter

#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom S4Vectors DataFrame
.add_values_to_colData <- function(x, values, name){
    # converts each value:name pair into a DataFrame
    values <- mapply(
        function(value, n){
            value <- DataFrame(value)
            colnames(value)[1L] <- n
            if(ncol(value) > 1L){
                i <- seq.int(2,ncol(value))
                colnames(value)[i] <- paste0(n,"_",colnames(value)[i])
            }
            value
        },
        values,
        name)

    values <- do.call(cbind, values)

    # check for duplicated values
    f <- colnames(colData(x)) %in% colnames(values)
    if(any(f)) {
        warning("The following values are already present in `colData` and ",
                "will be overwritten: '",
                paste(colnames(colData(x))[f], collapse = "', '"),
                "'. Consider using the 'name' argument(s) to specify alternative ",
                "names.",
                call. = FALSE)
    }
    # keep only unique values
    colData(x) <- cbind(colData(x)[!f], values)

    x
}

#' @importFrom stats acf pacf Box.test arima
# wrapper for time series measures
.getTsMeasures <- function(x, time.field, measure, ...) {
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
