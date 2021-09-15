#' Rank the sample information in a `SummarizedExperiment` Object
#'
#' The information stored in `colData` in \linkS4class{SummarizedExperiment}
#' object is ranked by `addRank` function and newly ranked field is added back
#' to \linkS4class{SummarizedExperiment} object.
#'
#' @param x \linkS4class{SummarizedExperiment} object
#' @param field A character value indicating the name of the matrix column that
#' is aimed to be ranked.
#' @param field_matrix A matrix that is aimed to be ranked. If the aimed matrix
#' string already exists in colData segment, this optional matrix doesn't need
#' to be added.
#' @param rank_field_name A character value to name the newly ranked `field`
#' @param ... Allow new parameters to be defined for this function.
#'
#' @return \linkS4class{RangedSummarizedExperiment} object with ranked field
#' added
#'
#' @examples
#' library(SummarizedExperiment)
#' data(airway, package="airway")
#' se <- airway
#'
#' abc <- addRank(se, field = "SampleName", rank_field_name = "SampleName_rank",
#'                                 na.last ="TRUE", ties.method = "first")
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom methods setGeneric
#'
#' @export
setGeneric("addRank", signature = "x",
        function(x, field, field_matrix, rank_field_name, ...)
            standardGeneric("addRank"))

setMethod("addRank", signature(x = "SummarizedExperiment"),
            function(x, field, field_matrix, rank_field_name, ...){
            col_data <- colData(x)
            mdat <- do.call(cbind,lapply(col_data, as.vector))

            if ( field %in% colnames(mdat) == TRUE){
                data <- mdat[, field]
                }  else {
                        data <- field_matrix
                        col_data[, field] <- data
            }

            ranked <- rank(data, ...)
            col_data[, rank_field_name] <- ranked
            colData(x) <- col_data
            return(x)
    })
