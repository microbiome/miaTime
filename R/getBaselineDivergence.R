#' Beta diversity between the baseline and later time steps 
#'
#' Calculates sample dissimilarity between the given baseline and other
#' time points, optionally within a group (subject, reaction chamber, or
#' similar). The corresponding
#' time difference is returned as well. The method operates on
#' `SummarizedExperiment` objects, and the results are stored in `colData`.
#'
#' @param x A
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#' @param group a single character value for specifying which grouping
#' factor is used (name of a `colData` field).
#' @param time_field a single character value, specifying the name of the
#' time series field in `colData`.
#' @param name_divergence a column vector showing beta diversity between samples
#' over n time intervals (default: \code{name_divergence = "time_divergence"})
#' @param name_timedifference field name for adding the time difference between
#' samples used to calculate beta diversity
#' (default: \code{name_timedifference = "time_difference"})
#' @param abund_values character indicating which assay values are used in
#' the dissimilarity estimation (default: \code{abund_values = "counts"})
#' @param FUN a \code{function} for dissimilarity calculation. The function must
#'   expect the input matrix as its first argument. With rows as samples 
#'   and columns as features. By default, \code{FUN} is
#'   \code{vegan::vegdist}.   
#' @param method a method that is used to calculate the distance. Method is
#'   passed to the function that is specified by \code{FUN}. By default,
#'   \code{method} is \code{"bray"}.
#'
#' @return a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' containing the sample dissimilarity and corresponding time difference between
#' samples (across n time steps), within each level of the grouping factor.
#'
#' @importFrom SEtools mergeSEs
#' @importFrom vegan vegdist
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
#'
#' @examples
#' #library(miaTime)
#' library(TreeSummarizedExperiment)
#'
#' data(hitchip1006)
#' tse <- mia::transformSamples(hitchip1006, method = "relabundance")
#'
#' # Subset to speed up example
#' tse <- tse[, colData(tse)$subject %in% c("900", "934", "843", "875")]
#'
#' tse2 <- getBaselineDivergence(tse,
#'                               group = "subject",
#'                               time_field = "time",
#'                               name_divergence = "divergence_from_baseline",
#'                               name_timedifference = "time_from_baseline",
#'                               abund_values="relabundance",
#'                               FUN = vegan::vegdist,
#'                               method="bray")
#'
#' @name getBaselineDivergence
#' @export
getBaselineDivergence <- function(x,
                            group=NULL,
                            time_field,
                            name_divergence = "divergence_from_baseline",
                            name_timedifference = "time_from_baseline",			    
                            abund_values = "counts",
			    FUN = vegan::vegdist,
			    method="bray"){

    # If group is not given, assume that all samples come from a single group
    if (is.null(group)) {
      spl <- split(seq_len(ncol(x)), rep(1, nrow(x)))
    } else {
      # Split SE into a list, by grouping
      if (is.factor(colData(x)[, group])) {
        colData(x)[, group] <- droplevels(colData(x)[, group])
      }
      spl <- split(seq_len(ncol(x)), colData(x)[, group])
    }

    # Apply the operation per subject
    x_list <- lapply(spl, function (s) {
        .calculate_divergence_from_baseline(x[,s], time_field, name_divergence, name_timedifference, abund_values, FUN, method)}
    )

    # Return the SE elements in a list
    if (length(x_list) > 1) {
        x2 <- mergeSEs(x_list)
    } else {
        x2 <- x_list[[1]]
    }

    # Just replace the colData for the original input
    colData(x) <- colData(x2)
    return(x)

}



# First define the function that calculates divergence for a given SE object
#' @importFrom mia estimateDivergence
.calculate_divergence_from_baseline <- function (x, time_field, name_divergence, name_timedifference, abund_values, FUN, method) {

    baseline_sample <- rownames(colData(x)[which.min(colData(x)[, time_field]),])

    # Add divergence from baseline
    base.sample <- as.vector(assay(x, abund_values)[, baseline_sample])
    x <- estimateDivergence(x, name = name_divergence, 
                     reference = base.sample,
                         FUN = FUN, method = method)

    # Add time difference from baseline
    colData(x)[, name_timedifference] <- colData(x)[, time_field] - colData(x)[baseline_sample, time_field]

    # Return
    return(x)

}


