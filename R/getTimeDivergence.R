#' Time difference between Nth time steps
#'
#' The time difference can be calculated between successive time steps, with
#' steps n >= 1.
#'
#' @param se A
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#' @param time_field column vector containing the time information
#' in `colData` field inside of \linkS4class{SummarizedExperiment} object
#' @param starting_position integer value that defines from which time step the
#' operation begins (default: \code{starting_position = 1})
#' @param time_interval integer value indicating the increment between the time
#' steps. It can be 1 or higher than 1.
#'
#' @return `colData` field in \linkS4class{SummarizedExperiment} object
#'
#' @importFrom SummarizedExperiment colData
#'
#' @examples
#' library(miaTime)
#' data(hitchip1006)
#' se <- hitchip1006
#'
#' seventeenTime <- getTimeDivergence(se, time_field = "time" , time_interval = 17)
#'
#'@export
getTimeDivergence <- function(se, time_field, starting_position = 1, time_interval){

    n <- 1:length(colData(se)[, time_field])

    if (!(time_interval %in% n)){
        stop("Time interval cannot be bigger than the length of time series")
    }

    vector <- colData(se)[, time_field]
    location <- seq(starting_position, length(vector), time_interval)

    x <- (vector)[location]

    time_diff <- c(NA, diff(x))

    timedivergence_n <- matrix(NA, ncol = 1, nrow = length(n))
    timedifference_n <- matrix(NA, ncol = 1, nrow = length(n))

    for (i in location){
        for (j in seq_len(length(location))){
            timedivergence_n[i] <- x[j]
            timedifference_n[i] <- time_diff[j]
        }
    }
    df <- as.data.frame(cbind(timedivergence_n, timedifference_n))
  return(df)
}
