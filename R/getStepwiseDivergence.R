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
#' @param time.col \code{Character scalar}. Specifies the name of the
#' time series field in `colData`.
#' @param time_interval \code{Integer scalar}. Indicates the increment between 
#' time steps. If you need to take every second, every third, or so, time step 
#' only, then increase this accordingly. (Default: \code{1})
#' @param name \code{Character scalar}. Shows beta diversity between 
#' samples. (Default: \code{"time_divergence"})
#' @param name.time \code{Character scalar}. Field name for adding the 
#' time difference between samples used to calculate beta diversity. 
#' (Default: \code{"time_difference"})
#' @param assay.type \code{Character scalar}. Specifies which assay values are 
#' used in the dissimilarity estimation. (Default: \code{"counts"})
#' @param method \code{Character scalar}. Used to calculate the distance. 
#' Method is passed to the function that is specified by \code{dis.fun}. 
#' (Default: \code{"bray"})
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
#'                               time.col = "time",
#'                               assay.type="relabundance",
#'                               dis.fun = vegan::vegdist,
#'                               method="bray")
NULL

#' @rdname addStepwiseDivergence
#' @export
#' 
#' @importFrom mia mergeSEs
#' @importFrom vegan vegdist
#' @importFrom mia addDivergence
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
        time.col,
        time_interval = 1,
        name = "time_divergence",
        name.time = "time_difference",
        assay.type = "counts",
        method="bray",
        n_dimred = NULL,
        dimred = NULL,
        ...){
        ############################# INPUT CHECK ##############################
        # name
        temp <- .check_input(
          name,
          list(NULL, "character scalar")
        )
        # name
        temp <- .check_input(
          name.time,
          list(NULL, "character scalar")
        )
        ########################### INPUT CHECK END ############################
        # Calculate values
        x <- .add_previous_sample(x, group, time.col, time_interval )
        res <- addDivergence(x, assay.type = assay.type, method = method, 
                             reference = "previous_sample", 
                             name = name, n_dimred = n_dimred, dimred = dimred, ...)
        col_data <- colData(res)
        colnames(col_data)[colnames(col_data) == "time_diff"] <- name.time
        colData(res) <- col_data 
        return(res)
    }
)

.add_previous_sample <- function(x, group, time, time_interval ){
  colData(x)$sample <- colnames(x)
  # For each group, get the same that has lowest time point
  df <- colData(x) %>% as.data.frame() %>%
    # Sort by subject and time
    arrange(.data[[group]], .data[[time]]) %>%
    group_by(.data[[group]]) %>%
    # Lag time by 1 (previous time point)
    mutate(previous_time = lag(.data[[time]], n = time_interval),  
           # Lag sample name by 1
           previous_sample = lag(sample, n = time_interval)) %>%  
    ungroup() |> DataFrame()
  rownames(df) <- df$sample
  df[["time_diff"]] <- df[[time]] - df[["previous_time"]]
  df <- df[ match(colnames(x), rownames(df)), ]
  colData(x) <- df
  return(x)
}
