#' Rank the sample information in a `SummarizedExperiment` Object
#'
#' The information stored in `colData` in \linkS4class{SummarizedExperiment}
#' object is ranked by `addRank` function and newly ranked field is added back
#' to \linkS4class{SummarizedExperiment} object.
#'
#' @param x \linkS4class{SummarizedExperiment} object
#' @param field A character string indicating the field aimed to be ranked
#' @param rank_field A single character value to name the newly ranked field
#' @param na.last Logical scalar: NA values are put last when \code{TRUE} ,first
#' when \code{FALSE}. If \code{'NA'}, they are removed; if \code{"keep"}
#' they are kept with rank NA.(default: \code{norm = TRUE})
#' @param ties.method a character string indicating the method used in ranking
#' @param ... Allow new parameters to be defined for this function.
#'
#' @return \linkS4class{SummarizedExperiment} object with ranked field added
#'
#' @examples
#' library(SummarizedExperiment)
#' data(airway, package="airway")
#' se <- airway
#'
#' addRank(se, field = "SampleName", rank_field = "SampleName_rank",
#'                                 na.last ="TRUE", ties.method = "first")
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom methods setGeneric
#'
#' @export
setGeneric("addRank", signature = "x",
        function(x, field, rank_field, na.last, ties.method, ...)
            standardGeneric("addRank"))

setMethod("addRank", signature(x = "SummarizedExperiment"),
            function(x, field, rank_field,
                    na.last = TRUE,
                    ties.method = c("average", "first", "last", "random",
                                                        "max", "min"), ...){
            col_data <- colData(x)
            mdat <- do.call(cbind,lapply(col_data, as.vector))
            data <- mdat[,field]
            ranked <- rank(data, na.last, ties.method)
            x$rank_field <- ranked
            return(x)
    })
