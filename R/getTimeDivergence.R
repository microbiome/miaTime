#' Beta diversity calculation between samples at the interval of n time steps
#'
#' The dissimilarity (beta distance) is calculated in steps of n for each
#' individual. "n" can be 1 or greater than 1. The time difference between
#' the time points used for beta diversity calculation is also given as output.
#'
#' @param se A
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object or
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' @param sample_field column vector containing the all individuals in
#' `colData` field
#' @param sample_id character indicating a individual inside of `sample_field`
#' @param time_field column vector containing the time information
#' in `colData` field
#' @param time_interval integer value indicating the increment between the time
#' steps. It can be 1 or higher than 1.
#' @param transposed logical scalar assigning samples to rows if they are in
#' columns (default: \code{transposed = FALSE})
#'
#' @return a list containing the beta diversity in matrix and a numeric vector
#' containing the time difference between the samples
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom vegan vegdist
#'
#' @examples
#' library(microbiomeDataSets)
#' tse <- GrieneisenTSData()
#'
#' BaboonDivergence <- getTimeDivergence(tse, sample_field = "baboon_id",
#'                                     sample_id = "Baboon_1" ,
#'                                     time_field = "collection_date",
#'                                     time_interval =2,
#'                                     transposed = FALSE )
#'
#'@export
getTimeDivergence <- function(se, sample_field, sample_id, time_field, time_interval, transposed = FALSE ){

    df <-as.data.frame(colData(se)[, c(sample_field, time_field)])

    #per sample
    sample_df <- df[which((df[, sample_field] == sample_id) == TRUE),]

    mat <- data.matrix(sample_df)

    #time in descending order
    mat[, time_field] <- mat[, time_field][order(mat[, time_field], decreasing = FALSE)]

    # Nth times are extracted
    location <-seq(0, length(mat[, time_field]), time_interval)

    extracted_time <- mat[, time_field][location]

    mat[, time_field] <- matrix(NA, ncol = 1, nrow = length(time))

    #extracted time places back
    mat[, time_field][location] <- extracted_time

    #related samples are chosen from assay
    x <- as.matrix(assay(tse)[,rownames(sample_df)[which(is.na(mat[, time_field]) == FALSE)]])

    if(!transposed){
      x <- t(x)
    }

    timedivergence_n <- as.matrix(vegan::vegdist(x), "bray")

    #timedivergence_n <-do.call(stats::dist, list(x))
    #timedivergence_n <- as.data.frame(as.matrix(time_divergence_n))

    time_diff <- diff(extracted_time)

    timedifference_n <- c(NA, time_diff)

    return(list(timedivergence_n = timedivergence_n,
                timedifference_n = timedifference_n))
}
