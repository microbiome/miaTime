#' Beta diversity calculation within individuals over time
#'
#' This method calculates sample dissimilarity between consecutive time steps
#' (step size n=1 by default), within a group (subject, reaction chamber,
#' or similar).
#' The corresponding time difference is returned as well. The method operates on
#' `SummarizedExperiment` objects, and the results are stored in `colData`.
#'
#' @param se A
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#' @param group a single character value for specifying which grouping
#' factor is used (name of a `colData` field).
#' @param time_field a single character value, specifying the name of the
#' time series field in `colData`.
#' @param time_interval integer value indicating the increment between time
#' steps (default: 1)
#' @param name_divergence a column vector showing beta diversity between samples
#' over n time intervals (default: \code{name_divergence = "time_divergence"})
#' @param name_timedifference field name for adding the time difference between
#' samples used to calculate beta diversity
#' (default: \code{name_timedifference = "time_difference"})
#' @param abund_values character indicating which assay values are used in
#' the dissimilarity estimation (default: \code{abund_values = "counts"})
#' @param distfun a function to calculate beta diversity.
#' (default: \code{distfun = vegan::vegdist})
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
#' @examples
#' #library(miaTime)
#' library(SummarizedExperiment)
#'
#' data(hitchip1006)
#' se <- mia::transformSamples(hitchip1006, method = "relabundance")
#'
#' # Subset to speed up example
#' se <- se[, colData(se)$subject %in% c("900", "934", "843", "875")]
#'
#' se2 <- getTimeDivergence(se, group = "subject",
#'                              time_interval = 1,
#'                              time_field = "time",
#'                              abund_values="relabundance",
#'                              distfun = vegan::vegdist)
#'
#' @name getTimeDivergence
#' @export
getTimeDivergence <- function(se,
                            group,
                            time_field,
                            time_interval=1,
                            name_divergence = "time_divergence",
                            name_timedifference = "time_difference",
                            abund_values = "counts",
			    distfun = vegan::vegdist){

    # Split SE into a list, by grouping
    spl <- split(colnames(se), colData(se)[, group])

    # Separate the groups with multiple time points
    spl_more <- spl[lapply(spl,length)>1]
    spl_one <- spl[lapply(spl,length) == 1]

    # Manipulate each subobject
    se_more_list <- lapply(seq_along(spl_more),
        function(i){.check_pairwise_dist(x = se[, spl_more[[i]]],
                                        distfun=distfun,
                                        time_interval,
                                        name_divergence = "time_divergence",
                                        name_timedifference = "time_difference",
                                        time_field,
                                        abund_values)})

            se_one_list <- lapply(seq_along(spl_one), function(i) {
                se[, spl_one[[i]]]}
            )

            for(i in seq_along(se_one_list)){
                colData(se_one_list[[i]])[, name_timedifference] <- NA
                colData(se_one_list[[i]])[, name_divergence] <- NA
            }

            # assign the names back to (T)SE objects
            names(se_more_list) <- names(spl_more)
            names(se_one_list) <- names(spl_one)

            # put lists together and put them in order
            whole_se <- do.call(c, list(se_one_list, se_more_list))
            whole_se <- whole_se[order(as.numeric(names(whole_se)))]

            # Merge the objects back into a single SE

            whole_se<-whole_se[!sapply(whole_se,is.null)]

            se_new <- mergeSEs(whole_se)

            if(class(whole_se[[1]]) == "TreeSummarizedExperiment"){
                colData(se) <- colData(se_new)
                return(se)
            }

            return(se_new)
}


.check_pairwise_dist <- function (x,
                                distfun,
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
                function (i) {distfun(mat[c(i, i-time_interval), ])})

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
