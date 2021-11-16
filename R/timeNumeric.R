#' Time field provided in character as "%H:%M:%S" or in Period class ("H M S")
#' can be converted to numeric values depending on the units.
#'
#' @param time_field : vector containing the time points
#' @param units : The unit of time points selected in `time_field`.
#' It can be hours, minutes or seconds.(default: \code{unit = "mins"})
#'
#' @examples
#' data(hitchip1006)
#' hitchipTime <- timeToDate(time_field = colData(hitchip1006)[,"time"],
#'                                                            unit = "mins")
#'
#' periodToNumeric <- timeNumeric(hitchipTime, unit = "mins")
#'
#' @export
timeNumeric <- function(time_field, units = "mins" ) {
    if (class(time_field) == "character"){

        hours <- substr(time_field, start = 1, stop = 2)
        minutes <- substr(time_field, start = 4, stop = 5)
        seconds <- substr(time_field, start = 7, stop = 8)

        num_h <- as.numeric(hours)
        num_m <- as.numeric(minutes)
        num_s <- as.numeric(seconds)

    } else if (class(time_field) == "Period"){

        num_h <- time_field$hour
        num_m <- time_field$minute
        num_s <- time_field$second
    }

    if (units == "hours"){
        time <- num_h + num_m/60 + num_s/3600

    } else if (units == "mins"){
        time <- num_h * 60 + num_m + num_s/60

    } else if (units == "secs"){
        time<- num_h * 3600 + num_m * 60 + num_s

    } else {
        warning('There is no valid time unit is given')
    }

    return(time)
}

