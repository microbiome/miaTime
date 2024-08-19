#' Beta diversity between the baseline and later time steps 
#'
#' Calculates sample dissimilarity between the given baseline and other
#' time points, optionally within a group (subject, reaction chamber, or
#' similar). The corresponding time difference is returned as well.
#' The method operates on `SummarizedExperiment` objects, and the results
#' are stored in `colData`.
#'
#' @inheritParams getStepwiseDivergence
#' @param baseline_sample \code{Character vector}. Specifies the baseline
#' sample(s) to be used. If the \code{group} argument is given, this must be a
#' named \code{vector}; one element per group.
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
#'
#' data(hitchip1006)
#' tse <- mia::transformAssay(hitchip1006, method = "relabundance")
#'
#' # Subset to speed up example
#' tse <- tse[, tse$subject %in% c("900", "934", "843", "875")]
#'
#' tse2 <- getBaselineDivergence(
#'     tse,
#'     group = "subject",
#'     time_field = "time",
#'     name_divergence = "divergence_from_baseline",
#'     name_timedifference = "time_from_baseline",
#'     assay.type="relabundance",
#'     FUN = vegan::vegdist,
#'     method="bray")
#'
#' tse2 <- getBaselineDivergence(
#'     tse,
#'     baseline_sample = "Sample-875",
#'     group = "subject",
#'     time_field = "time",
#'     name_divergence = "divergence_from_baseline",
#'     name_timedifference = "time_from_baseline",
#'     assay.type="relabundance",
#'     FUN = vegan::vegdist,
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

#' @rdname getPrevalence
#' @export
setMethod("addBaselineDivergence", signature = c(x = "SummarizedExperiment"),
    function(
        x,
        time_field,
        assay.type = "counts",
        group = NULL,
        name_divergence = "divergence_from_baseline",
        name_timedifference = "time_from_baseline",
        method = "bray",
        ...){
        ############################# INPUT CHECK ##############################
        # name_divergence
        temp <- .check_input(
            name_divergence,
            list(NULL, "character scalar")
        )
        # name_timedifference
        temp <- .check_input(
            name_timedifference,
            list(NULL, "character scalar")
        )
        ########################### INPUT CHECK END ############################
        # Calculate values
        res <- .get_baseline_divergence(
            x = x, group = group, time_field = time_field, assay.type = assay.type, method = method, ...)
        # Add values to colData
        x <- .add_values_to_colData(
            x, res, name = c(name_divergence, name_timedifference))
        return(x)
        
    }
)

.get_baseline_divergence <- function(
        x, group, baseline_sample, time_field, assay.type, method,
        altexp = NULL, baseline = NULL, ...){
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
    # time_field
    temp <- .check_input(
        time_field,
        list("character scalar"),
        supported_values = colnames(colData(x))
    )
    # Check that timepoints are numeric
    if( !is.numeric(x[[time_field]]) ){
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
    }
    # If not specified, for each group, get baseline sample. The baseline
    # sample is assumed to be a sample with lowest timepoint.
    if( is.null(baseline) ){
        baseline <- "baseline_sample"
        colData(x)[[baseline]] <- .get_baseline_sample(x, group, time_field)
    }
    # Check that baseline samples are correct
    .check_baseline_samples(x, baseline, group)
    ############################# INPUT CHECK END ##############################
    # Get a vector that shows which samples belong to which group 
    spl <- split(seq_len(ncol(x)), unfactor(colData(x)[[group]]))
    # Apply the operation per group; with group-specific baselines
    res <- lapply(names(spl), function(g){
        x_sub <- x[, spl[[g]]]
        res <- .calculate_divergence_from_baseline(
            x_sub, assay.type, method, time_field, baseline, add.ref = TRUE, ...)
        return(res)
        })
    res <- .wrangle_divergence_list(res, x)
    return(res)
}

.check_baseline_samples <- function(x, baseline, group){
    # Check that each group have only one baseline sample specified.
    baseline_samples <- split(colData(x)[[baseline]], unfactor(colData(x)[[group]]))
    correct <- lapply(baseline_samples, function(group){
        # Get unique
        group <- unique(group)
        # It must be a single index or character psacifying a column
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

.get_baseline_sample <- function(x, group, time){
    colData(x)$sample <- colnames(x)
    # For each group, get the sampe that has lowest time point
    baseline <- colData(x) %>% as.data.frame() %>%
        group_by(.data[[group]]) %>%
        mutate(rank = rank(.data[[time]], ties.method = "first")) %>%
        filter(rank == 1) %>%	
        select(.data[["sample"]], .data[[group]])
    # For each sample, assign corresponding baseline sample
    ind <- match(colData(x)[[group]], baseline[[group]])
    baseline <- baseline[ind, ]
    baseline <- baseline[["sample"]]
    return(baseline)
}

# First define the function that calculates divergence for a given SE object
#' @importFrom mia estimateDivergence
#' @importFrom methods is
.calculate_divergence_from_baseline <- function(
        x, assay.type, method, time_field, baseline,
        fun = FUN, FUN = vegan::vegdist, dimred = NULL, n_dimred = NULL, add.ref = TRUE, ...){
    # Get reference aka baseline sample
    ref_sample <- unique(x[[baseline]])
    
    reference <- x[, ref_sample]
    if( !add.ref ){
        not_ref <- colnames(x)[ !colnames(x) %in% ref_sample ]
        x <- x[, not_ref]
    }
    
    # Getting corresponding matrices, to calculate divergence 
    mat <- .get_mat_from_sce(x, assay.type, dimred, n_dimred)
    ref_mat <- .get_mat_from_sce(reference, assay.type, dimred, n_dimred)
    # transposing mat if taken from reducedDim. In reducedDim, samples are in
    # rows
    if( !is.null(dimred) ){
        mat <- t(mat)
        ref_mat <- t(ref_mat)
    }
    # Beta divergence from baseline info
    divergencevalues <- .calc_reference_dist(
        mat, as.vector(ref_mat), method, FUN = FUN, ...) ###################################### In mia, FUN --< stats::dist --> vegdist --> USE getDissimilarity????
    # Add time divergence from baseline info; note this has to be a list    
    timevalues <- colData(x)[[time_field]] - colData(reference)[[time_field]]
    names(divergencevalues) <- names(timevalues) <- colnames(x)
    
    res <- list(time = timevalues, divergence = divergencevalues)
    return(res)
}

.wrangle_divergence_list <- function(res, x){
    divergence <- lapply(res, function(values) values[["divergence"]])
    time <- lapply(res, function(values) values[["time"]])
    divergence <- unlist(divergence)
    time <- unlist(time)
    divergence <- divergence[ match(colnames(x), names(divergence))]
    time <- time[ match(colnames(x), names(time))]
    names(divergence) <- names(time) <- colnames(x)
    # Create a list of 2 elements. One element has all time differences, other
    # has all divergence values.
    res <- list(divergence, time)
    return(res)
}
