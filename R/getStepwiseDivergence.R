#' Beta diversity between consecutive time steps 
#'
#' Calculates sample dissimilarity between consecutive time points (t, t+i),
#' within a group (subject, reaction chamber, or similar). The corresponding
#' time difference is returned as well. The method operates on
#' `SummarizedExperiment` objects, and the results are stored in `colData`.
#'
#' @param x A
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#' @param group optional; a single character value for specifying the grouping
#' factor (name of a `colData` field). If given, the divergence is calculated
#' per group.  e.g. subject, chamber, group etc.).
#' @param time_field a single character value, specifying the name of the
#' time series field in `colData`.
#' @param time_interval integer value indicating the increment between time
#' steps (default: 1). If you need to take every second, every third, or so, time step only, then
#' increase this accordingly.
#' @param name_divergence a column vector showing beta diversity between samples
#' (default: \code{name_divergence = "time_divergence"})
#' @param name_timedifference field name for adding the time difference between
#' samples used to calculate beta diversity
#' (default: \code{name_timedifference = "time_difference"})
#' @param assay.type character indicating which assay values are used in
#' the dissimilarity estimation (default: \code{assay.type = "counts"}).
#' @param FUN a \code{function} for dissimilarity calculation. The function must
#'   expect the input matrix as its first argument. With rows as samples 
#'   and columns as features. By default, \code{FUN} is
#'   \code{vegan::vegdist}.
#' @param method a method that is used to calculate the distance. Method is
#'   passed to the function that is specified by \code{FUN}. By default,
#'   \code{method} is \code{"bray"}.
#' @param altexp String or integer scalar specifying the alternative experiment 
#' containing the input data.
#' @param dimred A string or integer scalar indicating the reduced dimension
#' result in `reducedDims` to use in the estimation.
#' @param n_dimred Integer scalar or vector specifying the dimensions to use if
#' \code{dimred} is specified.
#' @param ... Arguments to be passed
#'
#' @return a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' containing the sample dissimilarity and corresponding time difference between
#' samples (across n time steps), within each level of the grouping factor.
#'
#' @name addStepwiseDivergence
#' @export
#'
#' @examples
#' #library(miaTime)
#' library(TreeSummarizedExperiment)
#'
#' data(hitchip1006)
#' tse <- mia::transformCounts(hitchip1006, method = "relabundance")
#'
#' # Subset to speed up example
#' tse <- tse[, colData(tse)$subject %in% c("900", "934", "843", "875")]
#'
#' # Using vegdist for divergence calculation, one can pass
#' # the dissimilarity method from the vegan::vegdist options
#' # via the "method" argument
#' tse2 <- addStepwiseDivergence(tse, group = "subject",
#'                               time_interval = 1,
#'                               time_field = "time",
#'                               assay.type="relabundance",
#'                               FUN = vegan::vegdist,
#'                               method="bray")
NULL

#' @rdname addStepwiseDivergence
#' @export
#' 
#' @importFrom mia mergeSEs
#' @importFrom vegan vegdist
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom SingleCellExperiment altExp
#'
setGeneric("addStepwiseDivergence", signature = c("x"), function(x, ... )
  standardGeneric("addStepWi"))

#' @rdname addStepwiseDivergence
#' @export
setMethod("addStepwiseDivergence", signature = c(x = "ANY"),
    function(
        x,
        group=NULL,
        time_field,
        time_interval = 1,
        name_divergence = "time_divergence",
        name_timedifference = "time_difference",
        assay.type = "counts",
        method="bray",
        ...){
        ############################# INPUT CHECK ##############################
        # name_divergence
        temp <- .check_input(
          name_divergence,
          list(NULL, "character scalar")
        )
        # name_divergence
        temp <- .check_input(
          name_timedifference,
          list(NULL, "character scalar")
        )
        ########################### INPUT CHECK END ############################
        # Calculate values
        res <- .get_stepwise_divergence(
          x = x, group = group, time_field = time_field, time_interval = time_interval, assay.type = assay.type, method = method, ...)
        # Add values to colData
        x <- .add_values_to_colData(
          x, res, name = c(name_divergence, name_timedifference))
        return(x)
    }
)

.get_stepwise_divergence <- function(
    x,
    group=NULL,
    time_field,
    time_interval=1,
    name_divergence = "time_divergence",
    name_timedifference = "time_difference",
    assay.type = "counts",
    FUN = vegan::vegdist,
    method="bray",
    altexp = NULL,
    dimred = NULL,
    n_dimred = NULL,
    ...){
    ##########################################
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
    # time_interval
    temp <- .check_input(
      time_interval,
      list(NULL, "integer scalar")
    )
    ############################# INPUT CHECK END ##############################
    ### CAlculate values
    res <- lapply(colnames(x), function(sample){
        # Get previous time point
        res <- .calculate_divergence_from_prev_timepoint(
            x, sample, assay.type, method, time_field, time_interval, ...)
        return(res)
    })
    # Create a list of 2 elements. One element has all time differences, other
    # has all divergence values.
    res <- unlist(res, recursive = FALSE)
    return(res)
    
}

.get_previous_sample <- function(x, sample, time_field, time_interval){
    # Get values in this group
    x <- x[ , colData(x)[[group]] ==  colData(x[, sample])[[group]]]
    # Order
    x <- x[ , order(colData(x)[[time_field]])]
    # Get previous time point
    sample_ind <- which(colnames(x) == sample)
    prev_ind <- sample_ind - time_interval
    # If the value is possible return it
    res <- NA
    if( prev_ind > 0 ){
        res <- colnames(x)[prev_ind]
    }
    return(prev_ind)
}

.get_stepwise_divergence2 <- function(x,
                            group=NULL,
                            time_field,
                            time_interval=1,
                            name_divergence = "time_divergence",
                            name_timedifference = "time_difference",
                            assay.type = "counts",
                            FUN = vegan::vegdist,
                            method="bray",
                            altexp = NULL,
                            dimred = NULL,
                            n_dimred = NULL,
                            ...){

    # Store the original x
    xorig <- x
    
    # Use altExp if mentioned and available
    if (!is.null(altexp)) {
        .check_altExp_present(altexp, x)
        x <- altExp(x, altexp)
    }

    # Temporary sample ID
    x$tmp_sample_identifier_for_getStepwiseDivergence <- paste("SampleID", 1:ncol(x), sep="-")

    # If group is not given, assume that all samples come from a single group
    # TODO: switch to mia::splitOn    
    if (is.null(group)) {
      spl <- split(seq_len(ncol(x)), rep(1, nrow(x)))
    } else {
      # Split SE into a list, by grouping
      if (is.factor(colData(x)[, group])) {
        colData(x)[, group] <- droplevels(colData(x)[, group])
      }      
      spl <- split(seq_len(ncol(x)), colData(x)[, group])
    }

    # Separate the groups with multiple time points
    spl_more <- spl[lapply(spl,length) > 1]
    spl_one  <- spl[lapply(spl,length) == 1]

    # Manipulate each subobject
    x_more_list <- lapply(seq_along(spl_more),
        function(i){.check_pairwise_dist(x = x[, spl_more[[i]]],
                                        FUN=FUN,
                                        time_interval,
                                        name_divergence = name_divergence,
                                        name_timedifference = name_timedifference,
                                        time_field,
                                        assay.type,
                                        method,
                                        altexp,
                                        dimred,
                                        n_dimred)})

    x_one_list <- lapply(seq_along(spl_one), function(i) {
        x[, spl_one[[i]]]}
    )

    for(i in seq_along(x_one_list)){
        colData(x_one_list[[i]])[, name_timedifference] <- NA
        colData(x_one_list[[i]])[, name_divergence] <- NA
    }

    # assign the names back to (T)SE objects
    names(x_more_list) <- names(spl_more)
    names(x_one_list) <- names(spl_one)

    # put lists together and put them in order
    whole_x <- do.call(c, list(x_one_list, x_more_list))
    whole_x <- whole_x[order(as.numeric(names(whole_x)))]

    # Merge the objects back into a single X
    whole_x <- whole_x[!sapply(whole_x,is.null)]

    # Return the SE elements in a list
    if (length(whole_x) > 1) {
        x_new <- mergeSEs(whole_x)
    } else {
        x_new <- whole_x[[1]]
    }

    # Ensure that sample sorting matches between the input and output data
    inds <- match(x$tmp_sample_identifier_for_getStepwiseDivergence,
                  x_new$tmp_sample_identifier_for_getStepwiseDivergence)
    x_new <- x_new[, inds]

    # Add the new fields to colData
    # Just replace the colData for the original input
    # colData(xorig) <- colData(x_new)

    # Add beta divergence from baseline info; note this has to be a list
    timevalues <- list(colData(x_new)[, name_timedifference])
    divergencevalues <- list(colData(x_new)[, name_divergence])

    xorig <- .add_values_to_colData(xorig, timevalues, name_timedifference)
    xorig <- .add_values_to_colData(xorig, divergencevalues, name_divergence)    

    return(xorig)

}

.calculate_divergence_from_prev_timepoint <- function(
    x, sample, assay.type, method, time_field, time_interval,
    fun = FUN, FUN = vegan::vegdist, dimred = NULL, n_dimred = NULL, ...){
    # Get preiouv sample
    prev_sample <- .get_previous_sample(x, sample, time, time_interval)
    
    if( !is.na(prev_sample) ){
      
    }
    # If this sample has previous time step, calculate. Othwerwise five just NAs.
    # Calculate divergence between this timepints and prvious time point
    return(res)
}

.check_pairwise_dist <- function (x,
                                FUN,
                                time_interval,
                                name_divergence = "time_divergence",
                                name_timedifference = "time_difference",
                                time_field,
                                assay.type,
                                method,
                                altexp,
                                dimred,
                                n_dimred,
                                ...){

    mat <- .get_mat_from_sce(x, assay.type, dimred, n_dimred)
    ## transposing mat if taken from assay 
    if (is.null(dimred)) mat <- t(mat)
    
    time <- colData(x)[, time_field]

    ## Add new field to coldata
    colData(x)[, name_divergence]     <- rep(NA, nrow(mat))
    colData(x)[, name_timedifference] <- rep(NA, nrow(mat))

    if (nrow(mat) > time_interval) {

        ## beta diversity calculation
        n <- sapply(seq((time_interval+1), nrow(mat)), ## Do not use sapply
            function (i) {FUN(mat[c(i, i-time_interval), ], method=method, ...)}) ## MAYBE USE same method as in gtBaselineDivergence

        for(i in seq((time_interval+1), ncol(x))){
            colData(x)[, name_divergence][[i]] <- n[[i-time_interval]]
        }

        ## time difference calculation
        time <- sapply((time_interval+1):nrow(mat),
            function (i) {diff(colData(x)[c(i-time_interval, i), time_field])})

        for(i in seq((time_interval+1), nrow(colData(x)))){
            colData(x)[, name_timedifference][[i]] <- time[[i-time_interval]]
        }
    }
    return(x)

}
