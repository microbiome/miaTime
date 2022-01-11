#' Beta diversity calculation within individuals in the interval of n time steps
#'
#' The dissimilarity (beta diversity) is calculated in time steps of n, between
#' the samples of subject (individuals).Time interval can be 1 or greater than 1.
#' The time difference between n points is also calculated. Given input,
#' `SummarizedExperiment` or `TreeSummarizedExperiment` object, returns with
#' the calculations inserted in columns of `colData` field.
#'
#' @param se A
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object or
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' @param field a single character value for specifying which grouping
#' factor chosen from columns of `colData` field.
#' @param time_field a single character value for specifying time series
#' vector given in `colData` field.
#' @param time_interval integer value indicating the increment between the time
#' steps. It can be 1 or higher than 1. (default: \code{time_interval = 1})
#'
#' @return a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' containing the beta diversity and time difference between samples
#' in n time steps
#'
#' @importFrom SEtools mergeSEs
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
#'
#' @examples
#' library(miaTime)
#' data(hitchip1006)
#'
#' sub_hitchip <- hitchip1006[,800:1151]
#' hitchipTime_subset <- getTimeDivergence(sub_hitchip, field = "subject",
#'                                    time_interval = 1,
#'                                    time_field = "time")
#'
#' @rdname getTimeDivergence
#'
#' @export
check_pairwise_dist <- function (x,
                                 distfun,
                                 time_interval,
                                 new_field = "time_divergence",
                                 new_field2 = "time_difference",
                                 time_field,
                                 abund_values){

  mat <- t(assay(x, abund_values = "counts"))

  time <- colData(x)[, time_field]

  ## Add new field to coldata
  colData(x)[, new_field] <- rep(NA, nrow(mat))
  colData(x)[, new_field2] <- rep(NA, nrow(mat))

  if (nrow(mat) > time_interval) {

    ##beta diversity calculation
    n <- sapply((time_interval+1):nrow(mat),
                function (i) {distfun(mat[c(i, i-time_interval), ])})

    for(i in (time_interval+1):nrow(colData(x))){
      colData(x)[, new_field][[i]] <- n[[i-time_interval]]
    }

    ##time difference calculation
    time <- sapply((time_interval+1):nrow(mat),
                   function (i) {diff(colData(x)[c(i-time_interval, i), time_field])})

    for(i in (time_interval+1):nrow(colData(x))){
      colData(x)[, new_field2][[i]] <- time[[i-time_interval]]

    }
    return(x)
  }
}

#' @rdname getTimeDivergence
#' @name getTimeDivergence
#' @export
getTimeDivergence <- function(se, field, time_field, time_interval, new_field = "time_divergence",
                              new_field2 = "time_difference", abund_values = "counts" ){

# Split SE into a list, by subject
spl <- split(colnames(se), colData(se)[, field])

spl_more <- spl[lapply(spl,length)>1]

spl_one <- spl[lapply(spl,length) == 1]

# Manipulate each subobject
se_more_list <- lapply(seq_along(spl_more),
                  function(i){check_pairwise_dist(x = se[, spl_more[[i]]],
                                                   distfun=vegan::vegdist,
                                                   time_interval,
                                                   new_field = "time_divergence",
                                                   new_field2 = "time_difference",
                                                   time_field,
                                                   abund_values)})


se_one_list <- lapply(seq_along(spl_one), function(i) {
    se[, spl_one[[i]]]
})

# assign the names back to (T)SE objects
names(se_more_list) <- names(spl_more)
names(se_one_list) <- names(spl_one)

# put lists together and put them in order
whole_se <- do.call(c, list(se_one_list, se_more_list))
whole_se <- whole_se[order(as.numeric(names(whole_se)))]

# Merge the objects back into a single SE
se_new <- mergeSEs(whole_se)

return(se_new)

}

