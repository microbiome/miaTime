#' Numerically provided time points are converted to a Period class object in
#' "H M S" format based on unit. In addition to the time information,
#' date information can also be given.
#'
#' @param time_field : vector containing numeric time points
#' @param unit : character, indicating the unit of given times
#' (default: \code{unit = "mins"})
#' @param date : vector containing date information in "YYYY-MM-DD" format.
#' This information is optional.
#'
#' @examples
#' library(miaTime)
#' data(hitchip1006)
#'
#' ##without date information
#'
#' hitchipTime <- timeToDate(time_field = colData(hitchip1006)[,"time"], unit = "mins")
#'
#' ##with date
#'
#' DateColumn <- seq.Date(as.Date("2011-05-01"), by = "days", length.out = 1151)
#'
#' hitchipDate <- timeToDate(time_field = colData(hitchip1006)[,"time"] , unit = "mins",
#'                             date = DateColumn)
#'
#' ##This can be added back to \linkS4class{SummarizedExperiment} object that
#' #is used to store `hitchip1006` dataset.
#'
#' colData(hitchip1006)$clock <- hitchipTime
#'
#' @importFrom SummarizedExperiment colData<-
#' @importFrom lubridate hms
#'
#' @return a Period class object("H M S") or a list containing the date(s)
#' information provided and a Period class object containing time information
#'
#' @export
timeToDate<- function(time_field, unit = "mins", date){

    if (missing(date)){
      origin <- rep("1970-01-01", length(time_field))
    } else {
      origin <- as.Date(date)
    }
    if (unit == "mins"){
      hours <- 0
      minutes <- floor(time_field)
      seconds <- time_field - floor(time_field)
      seconds <- round(seconds, 2)*10
    } else if (unit == "hours"){
      hours <- floor(time_field)
      minutes <- time_field - floor(time_field)
      minutes <- round(minutes, 2)*10
      seconds <- 0
    } else if (unit == "secs"){
      hours <- 0
      minutes <- 0
      seconds <- floor(time_field)
    }

    if (unit == "mins"){
        ct <- as.POSIXct(60 * minutes, tz = "UTC", origin) + seconds
    } else if (unit == "hours"){
        ct <- as.POSIXct(3600 * hours, tz = "UTC", origin) + 60*minutes
    } else if (unit == "secs"){
        ct <- as.POSIXct(seconds, tz = "UTC", origin)
    } else {
        warning('The units must be hour(hours), minutes(mins), seconds(secs).
            The default value: mins is used')
    }

    clock <- format(ct, format="%H:%M:%S")
    period <- hms(clock)

    if (missing(date)){
      return(time = period)
    } else {
      return(list(time = period, date = origin))
    }

}
