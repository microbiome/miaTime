#' Beta diversity calculation between samples at the interval of n time steps
#'
#' The dissimilarity (beta diversity) is calculated in steps of n between the
#' samples of subject. "n" can be 1 or greater than 1. The time difference
#' between the n points is also calculated.
#'
#' @param se A
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object or
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' @param field column vector chosen from `colData` field to determine
#' the grouping factor in order to split `SummarizedExperiment` object.
#' @param time_field column vector chosen from `colData` field  indicating the
#' time
#' @param time_interval integer value indicating the increment between the time
#' steps. It can be 1 or higher than 1. (default: \code{time_interval = 1})
#'
#' @return a \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' containing the beta diversity and time difference between samples
#' in n time steps
#'
#' @examples
#' library(miaTime)
#' data(hitchip1006)
#'
#' hitchipTime <- getTimeDivergence(hitchip1006, field = "subject",
#'                                    time_interval = 1,
#'                                    time_field = "time")
#'
#'@export
getTimeDivergence <- function(se, field, time_field, time_interval){

# Define the names for the newly added fields
# for the time difference and community dissimilarity
new_field <- "time_divergence"
new_field2 <- "time_difference"

# Dissimilarity will be calculated between rows, as in stats::dist
check_pairwise_dist <- function (x, secondx, distfun, time_interval, new_field, new_field2, time_field) {

  if (nrow(colData(x)) == 1 | nrow(colData(secondx))== 1){
    mat <- t(assay(x))
    mat2 <- t(assay(secondx))
    mat_t <- rbind(mat, mat2)

    time <- colData(x)[, time_field]
    time2 <- colData(secondx)[, time_field]
    time_t <- rbind(time, time2)

    # Add new field to coldata
    colData(x)[, new_field] <- rep(NA, nrow(mat))

    if (nrow(mat_t) > time_interval) {
      d <- distfun(mat_t[c(2, 1), ])
      colData(x)[, new_field] <- d
      t <- diff(time_t)
      colData(x)[, new_field2] <- t
    }} else {
    colData(x)[, new_field] <- NA
    colData(x)[, new_field2] <- NA
  }
  return(x)
}

# Split SE into a list, by subject
spl <- split(colnames(se), colData(se)[, field])

# Manipulate each subobject
se_list <- lapply((time_interval+1):length(spl),
                  function(i){ check_pairwise_dist(x = se[, spl[[i - time_interval]]],
                                                   secondx = se[, spl[[i]]],
                                                   distfun=stats::dist,
                                                   time_interval,
                                                   new_field,
                                                   new_field2,
                                                   time_field)})

# Merge the objects back into a single SE
se <- do.call("cbind", se_list)
}

