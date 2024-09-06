#' @title Short term Changes in Abundance
#'
#' @description Calculates short term changes in abundance of taxa
#'              using temporal Abundance data.
#'  
#' @param x a \code{\link{SummarizedExperiment}} object.
#' 
#' @param assay.type \code{Character scalar}. Specifies the name of assay 
#'   used in calculation. (Default: \code{"counts"})
#'   
#' @param rarefy \code{Logical scalar}. Whether to rarefy counts.
#' (Default: \code{FALSE})
#' 
#' @param compositional \code{Logical scalar}. Whether to transform counts.
#' (Default: \code{FALSE})
#' 
#' @param depth \code{Integer scalar}. Specifies the depth used in rarefying. 
#' (Default: \code{min(assay(x, assay.type)}))
#' 
#' @param plot \code{Logical scalar}. Whether to plot short term change.
#' (Default: \code{FALSE})
#' 
#' @param ... additional arguments.
#' 
#' 
#' @return \code{dataframe} or \code{plot} with \code{short term change} 
#' calculations.
#' 
#' @details This approach is used by Wisnoski NI and colleagues
#'          \url{https://github.com/nwisnoski/ul-seedbank}. Their approach is based on
#'          the following calculation log(present abundance/past abundance).
#'          Also a compositional version using relative abundance similar to
#'          Brian WJi, Sheth R et al
#'          \url{https://www.nature.com/articles/s41564-020-0685-1} can be used.
#'          This approach is useful for identifying short term growth behaviors of taxa.
#'          
#' @name shortTermChange
#' 
#' 
#' @examples
#' 
#' # Load time series data
#' data(minimalgut)
#' tse <- minimalgut
#' 
#' short_time_labels <- c("74.5h", "173h", "438h", "434h", "390h")
#' 
#' # Subset samples by Time_lable and StudyIdentifier
#' tse <- tse[, !(colData(tse)$Time_label %in% short_time_labels)]
#' tse <- tse[, (colData(tse)$StudyIdentifier == "Bioreactor A")]
#' 
#' # Plot short term change in abundance
#' shortTermChange(tse, rarefy = TRUE, plot = TRUE) + ggtitle("Bioreactor A")


#' @rdname shortTermChange
#' @export
setGeneric("shortTermChange", signature = c("x"),
    function(x,assay.type = "counts", rarefy = FALSE, compositional = FALSE, 
        depth = min(assay(x, assay.type)), ...)
        standardGeneric("shortTermChange"))

#' @rdname shortTermChange
#' @export
#' @importFrom dplyr arrange as_tibble
#' @importFrom ggplot2 ggplot
#' @importFrom ggrepel geom_text_repel
#' @importFrom mia rarefyAssay
setMethod("shortTermChange", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = "counts", rarefy = FALSE, compositional = FALSE, 
        depth = min(assay(x, assay.type)), plot = FALSE, ...){
        ############################## Input check #############################
        # Check validity of object
        if(nrow(x) == 0L){
            stop("No data available in `x` ('x' has nrow(x) == 0L.)",
                 call. = FALSE)
        }
        # Check assay.type
        mia:::.check_assay_present(assay.type, x)
        if(!mia:::.is_a_bool(rarefy)){
            stop("'rarefy' must be TRUE or FALSE.", call. = FALSE)
        }
        if(!mia:::.is_a_bool(compositional)){
            stop("'compositional' must be TRUE or FALSE.", call. = FALSE)
        }
        # Ensure that the provided depth is valid
        if ( !is.null(depth) && depth > min(assay(x, assay.type)) ) {
            stop("Depth cannot be greater than the minimum number of counts in your data")
        
        }
        ############################ Input check end ###########################
        
        ############################ Data Preparation #########################
        # Initialize the filtered object based on rarefy and compositional arguments
        if (rarefy == TRUE && compositional == FALSE) {
            message("rarefy is set to TRUE, calculating short term change using counts")
            x <- rarefyAssay(x, assay.type = assay.type, depth = depth)
            assay.type <- "subsampled"
        } else if (rarefy == FALSE && compositional == FALSE) {
            message("rarefy is set to FALSE, compositional==FALSE, using raw counts")
            x <- x
        } else if (rarefy == FALSE && compositional == TRUE) {
            message("rarefy is set to FALSE, compositional==TRUE, using relative abundances")
            x <- transformAssay(x, method = "relabundance", assay.type = assay.type)
            assay.type <- "relabundance"
        } else if (rarefy == TRUE && compositional == TRUE) {
            stop("Both rarefy and compositional cannot be TRUE simultaneously")
        }
        ########################### Growth Metrics ############################
        grwt <- .calculate_growth_metrics(x, assay.type = assay.type, ...)
        #max growth
        maxgrs <- grwt %>%
            summarize(max.growth = max(growth_diff, na.rm = T))
        grs.all <- merge(grwt, maxgrs, by = "OTU")
        grs.all <- grs.all %>%
            mutate(ismax = ifelse(growth_diff == max.growth, T, F))

        grs.all$OTU <- gsub("_", " ", grs.all$OTU)
        
        grs.all$OTUabb <- toupper(abbreviate(grs.all$OTU,
                                             minlength = 3,
                                             method = "both.sides"
        ))
        grs.all$otu.time <- paste0(grs.all$OTUabb, " ", grs.all$time, "h")
        if(plot){
            grs.all <- ggplot(grs.all, aes(x = Time, group = OTU, col = OTU)) +
                geom_line(aes(y = growth_diff), alpha = 0.6, size = 1) +
                geom_point(
                    data = subset(grs.all, ismax == T),
                    aes(y = max.growth), alpha = 0.8, size = 3
                ) +
                geom_ribbon(aes(group = NULL, col = NULL, ymax = 0, ymin = -9),
                            fill = "#edf3f5", col = "#edf3f5", alpha = 0.5
                ) +
                theme(
                    legend.position = "right", legend.text = element_text(size = 9),
                    panel.background = element_rect(fill = "white"),
                    panel.grid.major = element_line(
                        size = 0.5, linetype = "solid",
                        colour = "#CCD1D1"
                    ),
                    panel.grid.minor = element_line(
                        size = 0.5, linetype = "solid",
                        colour = "#CCD1D1"
                    ),
                    legend.key = element_blank()
                ) +
                geom_text_repel(
                    data = subset(grs.all, ismax == T),
                    aes(y = max.growth, label = otu.time),
                    nudge_y = 0.2, box.padding = .5, max.iter = 10000,
                    size = 3, color = "black", segment.alpha = .5,
                    fontface = "italic", direction = "both"
                ) +
                geom_hline(yintercept = 0) +
                labs(
                    x = "Time (hr)",
                    y = expression(paste("Change in abundance ", " \U00B5 = ", abund[t + 1] / abund[t]))
                )
        }
        return(grs.all)
    
    }
)
        ########################################################################
# wrapper to calculate growth matrix
.calculate_growth_metrics <- function(x, assay.type, ...) {
    assay_data <- assay(x, assay.type)
    time <- colData(x)[, "Time.hr"]
    colnames(assay_data) <- time
    
    # Reshape data
    assay_data <- reshape2::melt(assay_data)
    assay_data <- as_tibble(assay_data) %>%
        arrange(Var2) %>%
        group_by(Var1) %>%
        mutate(
            time_lag = Var2 - lag(Var2), 
            growth_diff = value - lag(value),
            growth_rate = (value - lag(value)) / lag(value),
            var_abund = (value - lag(value)) / time_lag
        )
    
    colnames(assay_data)[colnames(assay_data) == "Var1"] <- "OTU"
    colnames(assay_data)[colnames(assay_data) == "Var2"] <- "Time"
    
    return(assay_data)
}