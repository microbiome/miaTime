#' Numerically provided time points are converted to a character object in the
#' format "%H:%M:%S".
#'
#' @param data : dataset with numeric time points in it
#' @param time_field : character, time points field in the dataset
#' @param unit : character, indicating the unit of given times
#' (default: \code{unit = "mins"})
#'
#' @examples
#' data(hitchip1006)
#' hitchipDate <- timeToDate(hitchip1006, time_field = "time", unit = "mins")
#'
#' ##This can be added back to \linkS4class{SummarizedExperiment} object that
#' #is used to store `hitchip1006` dataset.
#'
#' colData(hitchip1006)$clock <- hitchipDate
#'
#' @return a list object containing time in "%H:%M:%S" format.
#'
#' @export
timeToDate<- function(data, time_field, unit = "mins"){
    origin <- "1970-01-01"
    time <- colData(data)[,time_field]
    if (unit == "mins"){
      hours <- 0
      minutes <- floor(time)
      seconds <- time - floor(time)
      seconds <- round(seconds, 2)*10
    } else if (unit == "hours"){
      hours <- floor(time)
      minutes <- time - floor(time)
      minutes <- round(minutes, 2)*10
      seconds <- 0
    }

    if (unit == "mins"){
        ct <- as.POSIXct(60 * minutes, tz = "UTC", origin) + seconds
    } else if (unit == "hours"){
        ct <- as.POSIXct(3600 * hours, tz = "UTC", origin) + 60*minutes
    } else {
        warning('The units must be either hours or minutes.
            The default value: mins is used')
    }

    colData(data)$rank <- rank(colData(data)[,time_field])
    clock <- format(ct, format="%H:%M:%S")
    return(clock)
}
