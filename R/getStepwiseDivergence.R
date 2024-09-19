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
#' @param group \code{Character scalar}. Specifies the grouping
#' factor (name of a `colData` field). If given, the divergence is calculated
#' per group.  e.g. subject, chamber, group etc.). (Default: \code{NULL})
#' @param time_field \code{Character scalar}. Specifies the name of the
#' time series field in `colData`.
#' @param time_interval \code{Integer scalar}. Indicates the increment between 
#' time steps. If you need to take every second, every third, or so, time step 
#' only, then increase this accordingly. (Default: \code{1})
#' @param name_divergence \code{Character scalar}. Shows beta diversity between 
#' samples. (Default: \code{"time_divergence"})
#' @param name_timedifference \code{Character scalar}. Field name for adding the 
#' time difference between samples used to calculate beta diversity. 
#' (Default: \code{"time_difference"})
#' @param assay.type \code{Character scalar}. Specifies which assay values are 
#' used in the dissimilarity estimation. (Default: \code{"counts"})
#' @param FUN \code{Function} for dissimilarity calculation. The function must
#' expect the input matrix as its first argument. With rows as samples and 
#' columns as features. (Default: \code{vegan::vegdist})
#' @param method \code{Character scalar}. Used to calculate the distance. 
#' Method is passed to the function that is specified by \code{FUN}. 
#' (Default: \code{"bray"})
#' @param altexp \code{Character scalar} or \code{integer scalar}. Specifies the 
#' alternative experiment containing the input data. (Default: \code{NULL})
#' @param dimred \code{Character scalar} or \code{integer scalar}. indicates the 
#' reduced dimension result in `reducedDims` to use in the estimation. 
#' (Default: \code{NULL})
#' @param n_dimred \code{Integer vector}. Specifies the dimensions to use if
#' \code{dimred} is specified. (Default: \code{NULL})
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
#' tse <- mia::transformAssay(hitchip1006, method = "relabundance")
#'
#' # Subset to speed up example
#' tse <- tse[, colData(tse)$subject %in% c("900", "934", "843", "875")]
#'
#' # Using vegdist for divergence calculation, one can pass
#' # the dissimilarity method from the vegan::vegdist options
#' # via the "method" argument
#' tse <- addStepwiseDivergence(tse, group = "subject",
#'                               time_interval = 1,
#'                               time_field = "time",
#'                               assay.type="relabundance",
#'                               FUN = vegan::vegdist,
#'                               method="bray")
NULL

#' @rdname addStepwiseDivergence
#' @export
#' 
#' @importFrom mia mergeSEs getDivergence
#' @importFrom vegan vegdist
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom SingleCellExperiment altExp
#'
setGeneric("addStepwiseDivergence", signature = c("x"), function(x, ... )
  standardGeneric("addStepwiseDivergence"))

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
          x = x, group = group, time_field = time_field, 
          time_interval = time_interval, assay.type = assay.type, method = method, ...)
        # Add values to colData
        x <- .add_values_to_colData(
          x, res, name = c(name_divergence, name_timedifference))
        return(x)
    }
)

.get_stepwise_divergence <- function(
    x,
    group = NULL,
    time_field,
    time_interval = 1,
    name_divergence = "divergence",
    name_timedifference = "time_diff",
    assay.type = "counts",
    FUN = vegan::vegdist,
    method="bray",
    altexp = NULL,
    dimred = NULL,
    n_dimred = NULL,
    ...){
    ##########################################
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
    if( time_interval > ncol(x) ){
        stop("'time_interval' cannot be greater than the number of samples.", call. = FALSE)
    }
    # If TreeSE does not have column names, add
    if( is.null(colnames(x)) ){
      colnames(x) <- paste0("sample_", seq_len(ncol(x)))
    }
    # If group is not given, assume that all samples come from a single group
    if( !is.null(group) ){
        group <- "group"
        colData(x)[[group]] <- rep(1, nrow = nrow(x))
    }
    ############################# INPUT CHECK END ##############################
    
    # 1 Get previous sample for each sample.
    x <- .add_previous_sample(x, group, time_field, time_interval)
    res <- getDivergence(x, assay.type, method = method, 
              reference = "previous_sample", ...)
    res <- res <- list(res, x[["time_diff"]])
    return(res)
    
}

.add_previous_sample <- function(x, group, time, time_interval){
  colData(x)$sample <- colnames(x)
  # For each group, get the sampe that has lowest time point
  df <- colData(x) %>% as.data.frame() %>%
    # Sort by subject and time
    arrange(all_of(group), all_of(time)) %>%
    # Group by subject
    group_by(subject) %>%
    # Lag time by 1 (previous time point)
    mutate(previous_time = lag(time, n = time_interval),  
           # Lag sample name by 1
           previous_sample = lag(sample, n = time_interval)) %>%  
    ungroup() |> DataFrame()
  
  rownames(df) <- df$sample
  df[["time_diff"]] <- df[[time]] - df[["previous_time"]]
  df <- df[ match(colnames(x), rownames(df)), ]
  colData(x) <- df
  return(x)
}