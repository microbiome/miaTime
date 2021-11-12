#' Time field is provided in the format of "%H:%M:%S" can be converted to
#' numeric values depending on the units. Consecutive time difference of these
#' numeric values can be then calculated as well as the time difference based on
#' the certain time point with parameter `baseline`.
#'
#' @param time_field : character vector, contains the time points
#' @param baseline : integer, value to calculate the distance from certain time
#' point. (default: \code{baseline = 1})
#' @param units : The unit of time points selected in `time_field`.
#' It can be hours, minutes or seconds.(default: \code{unit = "mins"})
#'
#' @examples
#' data(hitchip1006)
#' hitchipDate <- timeToDate(hitchip1006, time_field = "time", unit = "mins")
#'
#' hitchipTime <- timeDifferences(hitchipDate, baseline = 1, units = "mins")
#'
#' ##Each of these fields in the list object can be applied to the dataset as
#' #new fields.
#'
#' colData(hitchip1006)$time_diff <- hitchipTime$Consecutive
#' colData(hitchip1006)$base_diff <- hitchipTime$Baseline
#'
#' @return a list object containing the numeric values of time, consecutive time
#' differences and the distance from the chosen time point to each time point.
#'
#' @export
timeDifferences <- function(time_field, baseline = 1, units = "mins"){
  hours <- substr(time_field, start = 1, stop = 2)
  minutes <- substr(time_field, start = 4, stop = 5)
  seconds <- substr(time_field, start = 7, stop = 8)

  num_h <- as.numeric(hours)
  num_m <- as.numeric(minutes)
  num_s <- as.numeric(seconds)

  if (units == "hours"){
    time <- num_h + num_m/60 + num_s/3600

  } else if (units == "mins"){
    time <- num_h * 60 + num_m + num_s/60

  } else if (units == "secs"){
    time<- num_h * 3600 + num_m * 60 + num_s

  } else {
    warning('There is no valid time unit input value is given')
  }

  diff_time <- diff(time)
  diff_time <- c(NA, diff_time)
  base <- time[baseline]
  i <- seq_len(length(time))
  ifelse(time[i], base_diff <- time[i] - base, base_diff )

  return(list(timeNumeric = time, Consecutive = diff_time,
              Baseline = base_diff))
}
