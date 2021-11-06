#'Numerically provided time points are converted to a POSIXct object in the
#'format "%Y-%m-%d %H:%M:%OS".
#'
#'@param data : dataset with time points in it
#'@param t : character time points field in the given dataset
#'@param origin : character indicating the start of the time points
#'(default: \code{origin = "1970-01-01"})
#'
#'@examples
#'data(hitchip1006)
#'force(hitchip1006)
#'HitchipDate <- TimePointsToDate(hitchip1006, "time", origin = "1970-01-01")
#'
#'@return a \linkS4class{TreeSummarizedExperiment} object containing added
#'date object, rank of the time points and time difference between the time
#'points in seconds.
#'
#'@export
TimePointsToDate<- function(data, t, origin = "1970-01-01"){
    time <- colData(data)[,t]
    seconds <- time - floor(time)

    minutes <- floor(time)
    i <- seq_len(length(seconds))
    ifelse(seconds[i] > 0.6, add_minutes <- seconds[i] %/% 0.6, add_minutes)
    minutes <- add_minutes + minutes
    seconds <- seconds - add_minutes*0.6

    seconds <- round(seconds, 2)*10
    ct <- as.POSIXct(60 * minutes, origin, tz = "UTC") + seconds
    colData(data)$Date <- ct
    colData(data)$rank <- rank(colData(data)[,t])
    z <- diff.POSIXt(ct)
    colData(data)$timeDifference<- c(NA,z)
    return(data)
}
