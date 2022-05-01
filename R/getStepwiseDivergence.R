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
#' @param ... Arguments to be passed
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
#' @alias getTimeDivergence
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
#' tse2 <- getStepwiseDivergence(tse, group = "subject",
#'                               time_interval = 1,
#'                               time_field = "time",
#'                               abund_values="relabundance",
#'                               FUN = vegan::vegdist)
#'
#' @name getStepwiseDivergence
#' @export
getStepwiseDivergence <- function(x,
                            group=NULL,
                            time_field,
                            time_interval=1,
                            name_divergence = "time_divergence",
                            name_timedifference = "time_difference",
                            abund_values = "counts",
			    FUN = vegan::vegdist, ...){


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

    # Separate the groups with multiple time points
    spl_more <- spl[lapply(spl,length) > 1]
    spl_one  <- spl[lapply(spl,length) == 1]

    # Manipulate each subobject
    x_more_list <- lapply(seq_along(spl_more),
        function(i){.check_pairwise_dist(x = x[, spl_more[[i]]],
                                        FUN=FUN,
                                        time_interval,
                                        name_divergence = "time_divergence",
                                        name_timedifference = "time_difference",
                                        time_field,
                                        abund_values)})

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

    # Just replace the colData for the original input
    colData(x) <- colData(x_new)
    
    return(x)

}




#' @rdname getStepwiseDivergence
#' @export
setGeneric("getTimeDivergence", signature = c("x"),
           function(x, ... )
             standardGeneric("getTimeDivergence"))

#' @rdname getStepwiseDivergence
#' @export
setMethod("getTimeDivergence",
          signature = c(x = "ANY"),
    function(x, ...){

        .Deprecated( msg = paste0("The name of the function 'getTimeDivergence' is",
                                  " changed to 'getStepwiseDivergence'. \nPlease use the new",
                                  " name instead.\n",
                                  "See help('Deprecated')") )

        getStepwiseDivergence(x, ...)
    }
)



.check_pairwise_dist <- function (x,
                                FUN,
                                time_interval,
                                name_divergence = "time_divergence",
                                name_timedifference = "time_difference",
                                time_field,
                                abund_values){

    mat <- t(assay(x, abund_values))

    time <- colData(x)[, time_field]

    ## Add new field to coldata
    colData(x)[, name_divergence]     <- rep(NA, nrow(mat))
    colData(x)[, name_timedifference] <- rep(NA, nrow(mat))

    if (nrow(mat) > time_interval) {

        ## beta diversity calculation
        n <- sapply((time_interval+1):nrow(mat),
                function (i) {FUN(mat[c(i, i-time_interval), ])})

        for(i in (time_interval+1):nrow(colData(x))){
                colData(x)[, name_divergence][[i]] <- n[[i-time_interval]]
        }

        ## time difference calculation
        time <- sapply((time_interval+1):nrow(mat),
            function (i) {diff(colData(x)[c(i-time_interval, i), time_field])})

        for(i in (time_interval+1):nrow(colData(x))){
            colData(x)[, name_timedifference][[i]] <- time[[i-time_interval]]
        }
    }
    return(x)

}
