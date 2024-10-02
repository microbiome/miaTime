#' Beta diversity between the baseline and later time steps 
#'
#' Calculates sample dissimilarity between the given baseline and other
#' time points, optionally within a group (subject, reaction chamber, or
#' similar). The corresponding time difference is returned as well.
#' The method operates on `SummarizedExperiment` objects, and the results
#' are stored in `colData`.
#'
#' @inheritParams getStepwiseDivergence
#' @param baseline_sample \code{Character vector}. Specifies the baseline sample(s) to be used. If the
#'   "group" argument is given, this must be a named vector; one element per group. (Default: \code{NULL})
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
#' The baseline sample/s always need to belong to the data object i.e. they can be merged into it before
#' applying this function. The reason is that they need to have comparable sample data, at least some time point
#' information for calculating time differences w.r.t. baseline sample.
#' 
#' The baseline time point is by default defined as the smallest time point (per group). Alternatively,
#' the user can provide the baseline vector, or a list of baseline vectors per group (named list per group).
#'
#' @importFrom mia mergeSEs
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom vegan vegdist
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom SingleCellExperiment altExp
#'
#' @examples
#' #library(miaTime)
#' library(TreeSummarizedExperiment)
#' library(dplyr)
#'
#' data(hitchip1006)
#' tse <- mia::transformCounts(hitchip1006, method = "relabundance")
#'
#' # Subset to speed up example
#' tse <- tse[, colData(tse)$subject %in% c("900", "934", "843", "875")]
#'
#' tse2 <- getBaselineDivergence(tse,
#'                               group = "subject",
#'                               time_field = "time",
#'                               name_divergence = "divergence_from_baseline",
#'                               name_timedifference = "time_from_baseline",
#'                               assay.type="relabundance",
#'                               FUN = vegan::vegdist,
#'                               method="bray")
#'
#' tse2 <- getBaselineDivergence(tse,
#'                               baseline_sample = "Sample-875",
#'                               group = "subject",
#'                               time_field = "time",
#'                               name_divergence = "divergence_from_baseline",
#'                               name_timedifference = "time_from_baseline",
#'                               assay.type="relabundance",
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
                            assay.type = "counts",
                            FUN = vegan::vegdist,
                            method="bray",
                            baseline_sample=NULL,
                            altexp = NULL,
                            dimred = NULL,
                            n_dimred = NULL,
                            ...){

    # Store the original data object
    xorig <- x
    
    # Use altExp if mentioned and available
    if (!is.null(altexp)) {
        .check_altExp_present(altexp, x)
        x <- altExp(x, altexp)
    }
    
    if (is.null(colnames(x))) {
        colnames(x) <- as.character(seq_len(ncol(x)))	
    }
    original.names <- colnames(x)

    # global vars
    is <- NULL    
    group_by <- NULL
    tmp_group_for_groupwise_splitting <- NULL
    time <- NULL
    filter <- NULL

    # Add time
    # colData(x)$time <- colData(x)[[time_field]]
    x <- .add_values_to_colData(x, list(colData(x)[[time_field]]), "time")

    # If group is not given, assume that all samples come from a single group    
    if (is.null(group)) {
        colData(x)$tmp_group_for_groupwise_splitting <- rep(1, nrow=nrow(x))
    } else if (is.character(group)) {
        colData(x)$tmp_group_for_groupwise_splitting <- as.character(colData(x)[[group]])
    } else {
        stop("The group argument in getBaselineDivergence should be NULL or a character i.e. name of a colData grouping field.")
    }

    # Split SE into a list, by grouping
    # TODO: switch to mia::splitOn
    spl <- split(seq_len(ncol(x)), colData(x)$tmp_group_for_groupwise_splitting)

    # Sample with the smallest time point within each subject
    # Use the smallest time point as the baseline
    if (is.null(baseline_sample)) {
        colData(x)$sample <- colnames(x)
        baseline <- colData(x) %>% as.data.frame() %>%
            group_by(tmp_group_for_groupwise_splitting) %>%
            mutate(rank = rank(time, ties.method="first")) %>%
	    filter(rank==1) %>%	
            select(sample, tmp_group_for_groupwise_splitting)
         baseline_sample <- baseline$sample
         names(baseline_sample) <- baseline$tmp_group_for_groupwise_splitting
	 nams <- names(baseline_sample)
      	 baseline_sample <- vapply(nams, function (g) {baseline_sample[[g]]}, "a")
	 names(baseline_sample) <- nams
    }

    # Then make sure that the baseline is an SE object
    if (is.character(baseline_sample)) {
         if (length(baseline_sample)==1) {
             baseline <- x[, baseline_sample]
         } else {
             if (is.null(names(baseline_sample))) {stop("Baseline sample has to be a named vector per group if it contains group-wise elements.")}
             # Just make sure that the given baseline samples are in the same order than the grouping variable
             baseline <- x[, baseline_sample[unique(colData(x)$tmp_group_for_groupwise_splitting)]]

	 }
    } else if (is(baseline_sample, "SummarizedExperiment")) {
        baseline <- baseline_sample
    } else {
        stop("Baseline sample not recognized in getBaselineDivergence. Should be NULL or a (named) character vector.")
    }

    # Apply the operation per group; with group-specific baselines
    if (ncol(baseline) == 1) {
        xli <- lapply(names(spl), function (g) {
            .calculate_divergence_from_baseline(x[,spl[[g]]], baseline,
	        time_field, name_divergence, name_timedifference, assay.type, FUN,
	        method, dimred, n_dimred, ...)})
    } else {
        xli <- lapply(names(spl), function (g) {
            .calculate_divergence_from_baseline(x[,spl[[g]]], baseline[, baseline_sample[[g]]],
	        time_field, name_divergence, name_timedifference, assay.type, FUN,
	        method, dimred, n_dimred, ...)})
    }

    # Return the elements in a list
    # FIXME: use SummarizedExperiment merge here or the new TreeSE merge thing
    if (length(xli) > 1) {
        x2 <- xli[[1]]
        for (i in seq(2, length(xli), 1)) {
            x2 <- TreeSummarizedExperiment::cbind(x2, xli[[i]])
	}
    } else {
        x2 <- xli[[1]]
    }

    # FIXME: reimplement the splitting so that we do not need intermediate variable like this
    colData(x2)$tmp_group_for_groupwise_splitting <- NULL

    # Return
    return(x2)

}


# First define the function that calculates divergence for a given SE object
#' @importFrom mia estimateDivergence
#' @importFrom methods is
.calculate_divergence_from_baseline <- function (x, baseline, time_field,
                                                 name_divergence, name_timedifference,
                                                 assay.type, FUN, method,
                                                 dimred, n_dimred) {

    # Global vars
    is <- NULL

    # If baseline is SE object then just ensure it has exactly one sample (well-defined baseline).
    # Otherwise, split the baseline from the data object.
    # Baseline is either an SE object with the same time field than x
    # or baseline specifies one sample from x
    if (is(baseline, "SummarizedExperiment")) {
        if (ncol(baseline)>1) {
            stop("If baseline is an SE object it should have a single sample.")
        } else {
            reference <- baseline
        }
    } else if (is.character(baseline) || is.numeric(baseline)) {
        reference <- x[, baseline]
    } else {
        stop("Baseline must be character or numeric vector specifying the SE sample; or it must be an SE object.")
    }

    # Getting corresponding matrices, to calculate divergence 
    mat <- .get_mat_from_sce(x, assay.type, dimred, n_dimred)
    ref_mat <- .get_mat_from_sce(reference, assay.type, dimred, n_dimred)
    
    # transposing mat if taken from reducedDim 
    if (!is.null(dimred)){
        mat <- t(mat)
        ref_mat <- t(ref_mat)
    }
    
    # Beta divergence from baseline info
    divergencevalues <- mia:::.calc_divergence(
        cbind(mat, ref_mat), colnames(ref_mat), FUN = FUN, method = method)
    divergencevalues <- divergencevalues[seq_len(ncol(mat)), "value"]

    # Add time divergence from baseline info; note this has to be a list    
    timevalues <- list(colData(x)[, time_field] - colData(reference)[, time_field])
    
    x <- .add_values_to_colData(x, timevalues, name_timedifference)
    x <- .add_values_to_colData(x, list(divergencevalues), name_divergence)    

    # Return
    return(x)

}


