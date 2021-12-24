#' Beta diversity calculation between samples at the interval of n time steps
#'
#' The dissimilarity (beta distance) is calculated in steps of n between the
#' samples of subject. "n" can be 1 or greater than 1. The time difference
#' between the n points is also calculated.
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
#' @return a matrix containing the beta diversity and time difference
#' between samples in n time steps
#'
#' @importFrom mia calculateDistance
#'
#' @examples
#' library(microbiomeDataSets)
#' tse <- GrieneisenTSData()
#'
#' BaboonDivergence <- getTimeDivergence(tse, sample_field = "baboon_id",
#'                                     sample_id = "Baboon_1" ,
#'                                     time_field = "collection_date",
#'                                     time_interval = 2,
#'                                     transposed = FALSE )
#'
#'@export
getTimeDivergence <- function(se, sample_field, sample_id, time_field, time_interval, transposed = FALSE ){

    df <-as.data.frame(colData(se)[, c(sample_field, time_field)])

    #per sample
    sample_df <- df[which((df[, sample_field] == sample_id) == TRUE),]

    mat <- data.matrix(sample_df)

    #time in increasing order
    mat[, time_field] <- mat[, time_field][order(mat[, time_field], decreasing = FALSE)]

     tt <- sapply(seq_len(nrow(mat)), FUN = function(i){
      k <- c(i, i + time_interval)
      i <- i + 1
      k
     })

    #total number of sample pairs
    total <- nrow(mat) - time_interval

    location <- tt[,seq_len(total)]

    #location matrix turned into list object
    list <- lapply(seq_len(ncol(location)), function(i) location[,i])

    samplename <-lapply(seq_len(length(list)), function(i){ mat[, sample_field][list[[i]]]})

    time <-lapply(seq_len(length(list)), function(i){ mat[, time_field][list[[i]]]})

    #related samples are chosen from assay
    list_t <- lapply(seq_len(length(time)), function(i){ as.matrix(assay(se)[,names(samplename[[i]])])})

    if(!transposed){
      list_t <- lapply(seq_len(length(list_t)), function(i){
        t(list_t[[i]])
        })
    }

    timedivergence_n <- vapply(seq_len(length(list_t)),
                               FUN = function(i) {calculateDistance(list_t[[i]])},
                               FUN.VALUE = 1)


    timedifference_n <- vapply(seq_len(length(timedivergence_n)),
                                     FUN = function(i) {diff(time[[i]])},
                                     FUN.VALUE = 1)

    return(cbind(timedivergence_n,timedifference_n))
}
