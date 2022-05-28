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
#' factor is used (name of a `colData` field). Optional.
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
#' @param baseline_sample Optional. The baseline sample to be used. A vector, or a named list if the
#' "group" argument is given. 
#'
#' @return a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' containing the sample dissimilarity and corresponding time difference between
#' samples (across n time steps), within each level of the grouping factor.
#'
#' @details
#' The group argument allows calculating divergence per group. Otherwise, this is done across all samples at once.
#' 
#' The baseline time point is by default defined as the smallest time point (per group). Alternatively,
#' the user can provide the baseline vector, or a list of baseline vectors per group (named list per group).
#'
#' @importFrom SEtools mergeSEs
#' @importFrom dplyr %>%
#' @importFrom vegan vegdist
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
#'
#' @examples
#' #library(miaTime)
#' library(TreeSummarizedExperiment)
#' library(dplyr)
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
#' tse2 <- getBaselineDivergence(tse,
#'                               baseline_sample = "Sample-875",
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
			    method="bray",
			    baseline_sample=NULL){

    # Store the original data object
    xorig <- x

    if (is.null(colnames(x))) {
        colnames(x) <- as.character(seq_len(ncol(x)))	
    }
    original.names <- colnames(x)

    # Add time
    colData(x)$time <- colData(x)[[time_field]]

    # If group is not given, assume that all samples come from a single group
    if (is.null(group)) {
        colData(x)$group <- rep(1, nrow=nrow(x))
    } else {
        colData(x)$group <- as.character(colData(x)[[group]])
    }

    # Split SE into a list, by grouping
    spl <- split(seq_len(ncol(x)), colData(x)$group)

    # Sample with the smallest time point within each subject
    # Use the smallest time point as the baseline
    if (is.null(baseline_sample)) {
        colData(x)$sample <- colnames(x)
        baseline <- colData(x) %>% as.data.frame() %>%
            group_by(group) %>%
            mutate(rank = rank(time, ties.method="first")) %>%
	    filter(rank==1) %>%	
            select(sample, group)
         baseline_sample <- baseline$sample
         names(baseline_sample) <- baseline$group
	 nams <- names(baseline_sample)
      	 baseline_sample <- lapply(nams, function (g) {x[, baseline_sample[[g]]]})
	 names(baseline_sample) <- nams
    } else if (is.character(baseline_sample)) {
         if (length(baseline_sample)==1) {
             baseline_sample <- list(x[, baseline_sample])
	 } else if (length(baseline_sample)>1) {
	     nams <- names(baseline_sample)
      	     baseline_sample <- lapply(nams, function (g) {x[, baseline_sample[[g]]]})
	     names(baseline_sample) <- nams
         }
    } else if (!class(baseline_sample)=="list") {
        baseline_sample <- list(baseline_sample)
    }

    # Apply the operation per group; with group-specific baselines
    if (length(baseline_sample) == 1) {
        x_list <- lapply(names(spl), function (g) {
            .calculate_divergence_from_baseline(x[,spl[[g]]], baseline_sample[[1]],
	        time_field, name_divergence, name_timedifference, abund_values, FUN, method)})
    } else {
        x_list <- lapply(names(spl), function (g) {
            .calculate_divergence_from_baseline(x[,spl[[g]]], baseline_sample[[g]],
	        time_field, name_divergence, name_timedifference, abund_values, FUN, method)})
    }

    # Return the SE elements in a list
    if (length(x_list) > 1) {
        x2 <- mergeSEs(x_list)
    } else {
        x2 <- x_list[[1]]
    }

    # Add the new fields (and only the new fields) into colData for the original input
    # Matching has to be done for cases that miss colnames. To be polished later.
    colData(xorig) <- cbind(colData(xorig), colData(x2)[match(original.names, colnames(x2)), c(name_timedifference, name_divergence)])

    xorig <- .add_values_to_colData(xorig, list(colData(x2)[match(original.names, colnames(x2)), name_timedifference]), name_timedifference)
    xorig <- .add_values_to_colData(xorig, list(colData(x2)[match(original.names, colnames(x2)), name_divergence]), name_divergence)    

    # Return
    return(xorig)

}


# First define the function that calculates divergence for a given SE object
#' @importFrom mia estimateDivergence
.calculate_divergence_from_baseline <- function (x, baseline, time_field, name_divergence, name_timedifference, abund_values, FUN, method, g) {

    # Reference vector
    reference <- as.vector(assay(baseline, abund_values))

    # Add beta divergence from baseline info
    x <- estimateDivergence(x, abund_values, name_divergence, reference, FUN, method)

    # Add time divergence from baseline info
    values <- list(colData(x)[, time_field] - colData(baseline)[, time_field])

    x <- .add_values_to_colData(x, values, name_timedifference)

    # Return
    return(x)

}


