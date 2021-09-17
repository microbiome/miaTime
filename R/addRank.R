#' Rank the sample information in a `TreeSummarizedExperiment` Object
#'
#' The information stored in `colData` in
#' \linkS4class{TreeSummarizedExperiment} object is ranked by `addRank`
#' function and newly ranked field is added back to
#' \linkS4class{TreeSummarizedExperiment} object.
#'
#' @param x \linkS4class{TreeSummarizedExperiment} object containing time series
#' @param field A character value indicating the name of the matrix column that
#' is aimed to be ranked.
#' @param rank_field_name A character value to name the newly ranked `field`
#' @param ... Allow new parameters to be defined for this function.
#'
#' @return \linkS4class{TreeSummarizedExperiment} object with ranked field
#' added
#'
#' @docType methods
#' @aliases addRank-TreeSummarizedExperiment
#' @aliases addRank,TreeSummarizedExperiment-method
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
#'
#' @examples
#' library(TreeSummarizedExperiment)
#' data(GlobalPatterns, package="mia")
#' se <- GlobalPatterns
#'
#' se1 <- addRank(se,field = "SampleType", rank_field_name = "SampleType_rank",
#'                 na.last = TRUE, ties.method = "first")
#'
#' #example 2: new matrix and its rank added
#'
#' se$newfield <- matrix(1:26, nrow = 26, ncol = 1)
#' se2 <- addRank(se, field = "newfield" ,
#'                 rank_field_name = "newfield_rank" , na.last = TRUE ,
#'                 ties.method = "first")
#'
#' @export
setGeneric("addRank", signature = "x",
        function(x, field, rank_field_name, ...)
            standardGeneric("addRank"))

setMethod("addRank", signature(x = "TreeSummarizedExperiment"),
            function(x, field, rank_field_name, ...){
            col_data <- colData(x)
            mdat <- do.call(cbind,lapply(col_data, as.vector))
            data <- mdat[, field]
            ranked <- rank(data, ...)
            col_data[, rank_field_name] <- ranked
            colData(x) <- col_data
            return(x)
    })
