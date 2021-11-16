#' Consecutive time difference of numeric values can be calculated
#' as well as the time difference based on the certain time point with
#' parameter `baseline`.
#'
#' @param time_field : vector containing time points
#' @param baseline : integer, value to calculate the distance from certain time
#' point. (default: \code{baseline = 1})
#'
#' @examples
#' data(hitchip1006)
#'
#' diffHitchip <- timeDifferences(time_field = colData(hitchip1006)[,"time"],
#'                                              baseline = 1)
#'
#' ##Each of these fields in the list object can be applied to the
#' #\linkS4class{SummarizedExperiment} object as new fields.
#'
#' colData(hitchip1006)$time_diff <- diffHitchip$Consecutive
#' colData(hitchip1006)$base_diff <- diffHitchip$Baseline
#'
#' @return a list object containing consecutive time differences and
#' the distance from selected time point to each time point.
#'
#' @export
timeDifferences <- function(time_field, baseline = 1){

    diff_time <- diff(time_field)
    diff_time <- c(NA, diff_time)
    base <- time_field[baseline]
    i <- seq_len(length(time_field))
    ifelse(time_field[i], base_diff <- time_field[i] - base, base_diff )

    return(list(Consecutive = diff_time,
                Baseline = base_diff))
}
