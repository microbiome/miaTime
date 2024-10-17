#' Beta diversity between the baseline and later time steps 
#'
#' Calculates sample dissimilarity between the given baseline and other
#' time points, optionally within a group (subject, reaction chamber, or
#' similar). The corresponding time difference is returned as well.
#' The method operates on `SummarizedExperiment` objects, and the results
#' are stored in `colData`.
#'
#' @inheritParams addStepwiseDivergence
#' @param dis.fun \code{Function} for dissimilarity calculation. The function must
#' expect the input matrix as its first argument. With rows as samples and 
#' columns as features. (Default: \code{vegan::vegdist})
#' @param baseline_sample \code{Character vector}. Specifies the baseline
#' sample(s) to be used. If the \code{group} argument is given, this must be a
#' named \code{vector}; one element per group.
#' @param ... optional arguments
#'
#' @return a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' containing the sample dissimilarity and corresponding time difference between
#' samples (across n time steps), within each level of the grouping factor.
#'
#' @details
#' The group argument allows calculating divergence per group. Otherwise, this
#' is done across all samples at once.
#'
#' The baseline sample/s always need to belong to the data object i.e. they
#' can be merged into it before
#' applying this function. The reason is that they need to have comparable
#' sample data, at least some time point
#' information for calculating time differences w.r.t. baseline sample.
#' 
#' The baseline time point is by default defined as the smallest time point
#' (per group). Alternatively,
#' the user can provide the baseline vector, or a list of baseline vectors per
#' group (named list per group).
#'
#' @examples
#' library(miaTime)
#' library(mia)
#'
#' data(hitchip1006)
#' tse <- transformAssay(hitchip1006, method = "relabundance")
#'
#' # Subset to speed up example
#' tse <- tse[, tse$subject %in% c("900", "934", "843", "875")]
#'
#' tse2 <- addBaselineDivergence(
#'     tse,
#'     group = "subject",
#'     time.col = "time",
#'     name = "divergence_from_baseline",
#'     name.time = "time_from_baseline",
#'     assay.type="relabundance",
#'     dis.fun = vegan::vegdist,
#'     method="bray")
#'
#' tse2 <- addBaselineDivergence(
#'     tse,
#'     baseline_sample = "Sample-875",
#'     group = "subject",
#'     time.col = "time",
#'     name = "divergence_from_baseline",
#'     name.time = "time_from_baseline",
#'     assay.type="relabundance",
#'     dis.fun = vegan::vegdist,
#'     method="bray")
#'
#' @name addBaselineDivergence
#' @export
#' 
NULL

#' @rdname addBaselineDivergence
#' @export
#' 
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
setGeneric("addBaselineDivergence", signature = "x", function(x, ...)
    standardGeneric("addBaselineDivergence"))

#' @rdname addBaselineDivergence
#' @export
setMethod("addBaselineDivergence", signature = c(x = "SummarizedExperiment"),
    function(
        x,
        time.col,
        assay.type = "counts",
        group = NULL,
        name = "divergence",
        name.time = "time_diff",
        method = "bray",
        dimred = NULL,
        ndimred = NULL,
        dis.fun = vegan::vegdist,
        baseline_sample = NULL,
        ...){
        ############################# INPUT CHECK ##############################
        # name
        temp <- .check_input(
            name,
            list(NULL, "character scalar")
        )
        # name.time
        temp <- .check_input(
            name.time,
            list(NULL, "character scalar")
        )
        ########################### INPUT CHECK END ############################
        # Calculate values
        x <- .get_baseline_divergence( x = x, group = group, 
                                       time.col = time.col, 
                                       assay.type = assay.type, 
                                       method = method, 
                                       name = name,
                                       name.time = name.time, 
                                       dimred = dimred, ndimred = ndimred, 
                                       baseline_sample = baseline_sample, ...)
        
        return(x)
        
    }
)

.get_baseline_divergence <- function(
        x, group, baseline_sample = NULL, 
        time.col, assay.type, method,
        altexp = NULL, baseline = NULL, 
        dimred = NULL, ndimred = NULL, 
        dis.fun = vegan::vegdist, 
        name.time = "time_diff", 
        name = "divergence", ...){
    ############################### INPUT CHECK ################################
    # If TreeSE does not have column names, add
    if( is.null(colnames(x)) ){
        colnames(x) <- as.character(seq_len(ncol(x)))	
    }
    # Use altExp if mentioned and available
    if( !is.null(altexp) ){
        .check_altExp_present(altexp, x)
        x <- altExp(x, altexp)
    }
    # assay.type
    .check_assay_present(assay.type, x)
    # time.col
    temp <- .check_input(
        time.col,
        list("character scalar"),
        supported_values = colnames(colData(x))
    )
    # Check that timepoints are numeric
    if( !is.numeric(x[[time.col]]) ){
        stop("Timepoints must be numeric.", call. = FALSE)
    }
    # group
    temp <- .check_input(
        group,
        list(NULL, "character scalar"),
        supported_values = colnames(colData(x))
    )
    # baseline
    temp <- .check_input(
        baseline,
        list(NULL, "character scalar"),
        supported_values = colnames(colData(x))
    )
    # If group is not given, assume that all samples come from a single group
    if( is.null(group) ){
        group <- "group"
        colData(x)[[group]] <- rep(1, nrow = nrow(x))
    } else if (is.character(group)) {
        colData(x)$group <- as.character(colData(x)[[group]])
    } else {
        stop("The group argument in getBaselineDivergence should be 
             NULL or a character i.e. name of a colData grouping field.")
    }
    # If not specified, for each group, get baseline sample. The baseline
    # sample is assumed to be a sample with lowest timepoint.
    # Sample with the smallest time point within each subject
    # Use the smallest time point as the baseline
    if (is.null(baseline_sample)) {
        colData(x)$sample <- colnames(x)
        baseline <- colData(x) %>% as.data.frame() %>%
            group_by(group) %>%
            mutate(rank = rank(time.col, ties.method="first")) %>%
            filter(rank==1) %>%	
            select(sample, group)
        baseline_sample <- baseline$sample
        names(baseline_sample) <- baseline$group
        nams <- names(baseline_sample)
        baseline_sample <- vapply(nams, function (g) {baseline_sample[[g]]}, "a")
        names(baseline_sample) <- nams
    }
    
    # Then make sure that the baseline is an SE object
    if (is.character(baseline_sample)) {
        if (length(baseline_sample)==1) {
            baseline <- x[, baseline_sample]
        } else {
            if (is.null(names(baseline_sample))) {stop("Baseline sample has to 
                be a named vector per group if it contains group-wise elements.")}
    # Just make sure that the given baseline samples are in the same order than 
    # the grouping variable
            baseline <- x[, baseline_sample[unique(colData(x)$group)]]
            
        }
    } else if (is(baseline_sample, "SummarizedExperiment")) {
        baseline <- baseline_sample
    } else {
        stop("Baseline sample not recognized in getBaselineDivergence. 
             Should be NULL or a (named) character vector.")
    }
    # Check that baseline samples are correct
    # .check_baseline_samples(x, baseline, group)
    ############################# INPUT CHECK END ##############################
    # Apply the operation per group; with group-specific baselines
    spl <- split(seq_len(ncol(x)), colData(x)$group)
    
    if (ncol(baseline) == 1) {
        xli <- lapply(names(spl), function (g) {
            .calculate_divergence_from_baseline(x[,spl[[g]]], baseline,
                                                time.col, name, 
                                                name.time, 
                                                assay.type, dis.fun,
                                                method, dimred, ndimred, ...)})
    } else {
        xli <- lapply(names(spl), function (g) {
            .calculate_divergence_from_baseline(x[,spl[[g]]], 
                                                baseline[, baseline_sample[[g]]],
                                                time.col, 
                                                name, 
                                                name.time, 
                                                assay.type, dis.fun,
                                                method, dimred, ndimred, ...)})
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
    
    # FIXME: reimplement the splitting so that we do not need intermediate 
    # variable like this
    colData(x2)$group <- NULL
    
    # Return
    return(x2)
}

.check_baseline_samples <- function(x, baseline, group){
    # Check that each group have only one baseline sample specified.
    baseline_samples <- split(colData(x)[[baseline]], unfactor(colData(x)[[group]]))
    correct <- lapply(baseline_samples, function(group){
        # Get unique
        group <- unique(group)
        # It must be a single index or character specifying a column
        res <- length(group) == 1 && (
            (is.integer(group) && group >= 1 && group <= ncol(x)) ||
                (is.character(group) && group %in% colnames(x)) )
        return(res)
    })
    correct <- unlist(correct)
    if( !all(correct) ){
        stop(
            "Each group must have only one baseline sample specified. ",
            "Moreover the 'baseline' must specify an index or name that ",
            "points to a column.", call. = FALSE)
    }
    return(NULL)
}

# First define the function that calculates divergence for a given SE object
#' @importFrom methods is
.calculate_divergence_from_baseline <- function (x, baseline, time.col,
                                                 name, 
                                                 name.time,
                                                 assay.type, dis.fun, method,
                                                 dimred, ndimred) {
    
    # If baseline is SE object then just ensure it has exactly one sample 
    # (well-defined baseline).
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
        stop("Baseline must be character or numeric vector specifying the 
             SE sample; or it must be an SE object.")
    }
    
    # Getting corresponding matrices, to calculate divergence 
    mat <- .get_mat_from_sce(x, assay.type, dimred, ndimred)
    ref_mat <- .get_mat_from_sce(reference, assay.type, dimred, ndimred)
    
    # transposing mat if taken from reducedDim 
    if (!is.null(dimred)){
        mat <- t(mat)
        ref_mat <- t(ref_mat)
    }
    
    # Beta divergence from baseline info
    divergencevalues <- mia:::.calc_divergence(
        cbind(mat, ref_mat), colnames(ref_mat), dis.fun = dis.fun, method = method)
    divergencevalues <- divergencevalues[seq_len(ncol(mat)), "value"]
    
    # Add time divergence from baseline info; note this has to be a list    
    timevalues <- list(colData(x)[, time.col] - colData(reference)[, time.col])
    
    x <- .add_values_to_colData(x, timevalues, name.time)
    x <- .add_values_to_colData(x, list(divergencevalues), name)    
    
    # Return
    return(x)
    
}
