setGeneric("addRank", signature = c("x"),
           function(x,
                    method = c("average", "first", "last", "random", "max", "min"))
             standardGeneric("addRank"))

setMethod("addRank", signature = c(x = "SummarizedExperiment"),
          function(x,
                   method = c("average", "first", "last", "random", "max", "min")){
          col_data <- colData(x)
          if (!is.vector(col_data)){
            col_data <- as.vector(col_data)
          }

          ranked <- rank(col_data, na.last = TRUE , method)
          return(ranked)
          }
)
