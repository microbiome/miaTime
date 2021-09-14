#' Rank the sample information in a `SummarizedExperiment` Object
#'
#' The information stored in `colData` in \linkS4class{SummarizedExperiment}
#' object can be ranked by `addRank` function.
#'
#' @param x \linkS4class{SummarizedExperiment} object
#' @param ... Allow new parameters to be defined for this function.
#' @param na.last  NA values are put last when \code{TRUE} or put first when
#' \code{FALSE} or they are kept with rank NA when \code{"keep"}
#' @param ties.method a character string indicating the method used in ranking
#'
#' @return a list shows the ranked components of
#' \linkS4class{SummarizedExperiment} object
#'
#' @examples
#' library(SummarizedExperiment)
#' data(airway, package="airway")
#' se <- airway
#'
#' addRank(se, na.last ="TRUE", ties.method = "last")
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @export
setGeneric("addRank", signature = "x",
        function(x,...)
            standardGeneric("addRank"))

.is_integer <- function(t){
    i <- seq_along(t)
    for (i in t){
        if(is.integer(t[[i]]) == FALSE){
        t[-i]
        }
    }
}

setMethod("addRank", signature(x = "SummarizedExperiment"),
            function(x,
                    na.last = "TRUE",
                    ties.method = c("average", "first", "last", "random", "max",
                        "min"), ...){
            col_data <- colData(x)
            list <- col_data@listData
            integerlist <- .is_integer(list) #lapply(list, .is_integer)
            ranked <- rank(integerlist, na.last, ties.method)
            return(ranked)
    })
