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
#' @importFrom mia mergeSEs
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
    # 2 Calculate dissimilarity matrix
    mat <- assay(x, assay.type)
    mat <- t(mat)
    res <- vegdist(mat, method = "bray")
    res <- as.matrix(res)
    # 3 Assign divergence based on dissimilarity matrix and previous sample information.
    mapping <- data.frame(sample = x$sample, prev_sample = colData(x)[["previous_sample"]])
    res <- mapping %>%
      rowwise() %>%
      mutate(divergence = get_divergence(sample, prev_sample, res)) %>%
      ungroup()
    res <- res[ match(res$sample, colnames(x)), ]
    x[["divergence"]] <- res[["divergence"]]
    res <- list(x[["previous_time"]], x[["divergence"]])
    return(res)
    
}

get_divergence <- function(current_sample, previous_sample, dissim_matrix) {
  if (!is.na(previous_sample) && 
      current_sample %in% rownames(dissim_matrix) && 
      previous_sample %in% colnames(dissim_matrix)) {
    return(dissim_matrix[current_sample, previous_sample])
  } else {
    return(NA)
  }
}

.calculate_divergence_based_on_reference <- function(
    x, assay.type, method, time_field, baseline,
    fun = FUN, FUN = vegan::vegdist, dimred = NULL, n_dimred = NULL, add.ref = TRUE, ...){
  # Get reference aka baseline sample
  prev_samples <- colData(x)[["previous_sample"]]
  prev_samples <- prev_samples[ !is.na(prev_samples) ]
  prev_samples <- x[ , x$sample %in% prev_samples ]
  # Getting corresponding matrices, to calculate divergence 
  mat <- .get_mat_from_sce(x, assay.type, dimred, n_dimred)
  ref_mat <- .get_mat_from_sce(prev_samples, assay.type, dimred, n_dimred)
  # transposing mat if taken from reducedDim. In reducedDim, samples are in
  # rows
  if( !is.null(dimred) ){
    mat <- t(mat)
    ref_mat <- t(ref_mat)
  }
  # Beta divergence from baseline info
  divergencevalues <- .calc_reference_dist(
    mat, ref_mat, method, FUN = FUN, ...) ###################################### In mia, FUN --< stats::dist --> vegdist --> USE getDissimilarity????
  # Add time divergence from baseline info; note this has to be a list    
  timevalues <- colData(x)[[time_field]] - colData(reference)[[time_field]]
  names(divergencevalues) <- names(timevalues) <- colnames(x)
  
  res <- list(time = timevalues, divergence = divergencevalues)
  return(res)
}

.calculate_divergence_from_prev_timepoint <- function(
        x, sample, assay.type, method, time_field, time_interval, group,
        fun = FUN, FUN = vegan::vegdist, dimred = NULL, n_dimred = NULL, ...){
    # Get preiouv sample
    prev_sample <- .get_previous_sample(x, sample, time_field, time_interval, group)
    # res <- list(time = NA, divergence = NA)
    # if( !is.na(prev_sample) ){
    #   x <- x[, c(sample, prev_sample)]
    #   ref_name <- "reference"
    #   colData(x)[[ref_name]] <- prev_sample
    #   res <- .calculate_divergence_from_baseline(x, assay.type, method, time_field, ref_name, add.ref = FALSE)
    # }
    # # If this sample has previous time step, calculate. Othwerwise, give just NAs.
    # # Calculate divergence between this timepints and prvious time point
    return(prev_sample)
}

.get_previous_sample <- function(x, sample, time_field, time_interval, group){
    # Get values in this group
    x <- x[ , colData(x)[[group]] ==  colData(x[, sample])[[group]]]
    # Order
    x <- x[ , order(colData(x)[[time_field]])]
    # Get previous time point
    sample_ind <- which(colnames(x) == sample)
    sample_ind <- sample_ind[[1]]
    prev_ind <- sample_ind - time_interval
    # If the index is possible, return sample name. Otherwise, return NA.
    res <- NA
    if( prev_ind > 0 ){
        res <- colnames(x)[prev_ind]
    }
    return(res)
}

.add_previous_sample <- function(x, group, time, time_interval){
  colData(x)$sample <- colnames(x)
  # For each group, get the sampe that has lowest time point
  df <- colData(x) %>% as.data.frame() %>%
    arrange(.data[[group]], .data[[time]]) %>%                 # Sort by subject and time
    group_by(subject) %>%                      # Group by subject
    mutate(previous_time = lag(time, n = time_interval),   # Lag time by 1 (previous time point)
           previous_sample = lag(sample, n = time_interval)) %>%  # Lag sample name by 1
    ungroup() |> DataFrame()
  rownames(df) <- colnames(x)
  colData(x) <- df
  return(x)
}
