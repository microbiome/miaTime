#' @rdname shortTermChange
#' @export
setGeneric("shortTermChange", signature = c("x"),
    function(x,assay.type = "counts", rarefy = FALSE, compositional = FALSE, 
        depth = NULL)
        standardGeneric("shortTermChange"))

#' @rdname shortTermChange
#' @export
setMethod("shortTermChange", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = "counts", rarefy = FALSE, compositional = FALSE, 
        depth = NULL){
        ############################## Input check #############################
        # Check validity of object
        if(nrow(x) == 0L){
            stop("No data available in `x` ('x' has nrow(x) == 0L.)",
                 call. = FALSE)
        }
        # Check assay.type
        .check_assay_present(assay.type, x)
        if(!.is_a_bool(rarefy)){
            stop("'rarefy' must be TRUE or FALSE.", call. = FALSE)
        }
        if(!.is_a_bool(compositional)){
            stop("'compositional' must be TRUE or FALSE.", call. = FALSE)
        }
        # Ensure that the provided depth is valid
        if (!is.null(depth)) {
            if (depth > min(assay(x, assay.type))) {
                stop("Depth cannot be greater than the minimum number of counts in your data")
            }
        } else {
            stop("Depth cannot be NULL")
        }
        
        ############################ Input check end ###########################
        
        ############################ Data Preparation #########################
        # Initialize the filtered object based on rarefy and compositional arguments
        if (rarefy == TRUE && compositional == FALSE) {
            message("rarefy is set to TRUE, calculating short term change using counts")
            x <- rarefyAssay(x, assay.type = assay.type, depth = depth)
        } else if (rarefy == FALSE && compositional == FALSE) {
            message("rarefy is set to FALSE, compositional==FALSE, using raw counts")
            x <- x
        } else if (rarefy == FALSE && compositional == TRUE) {
            message("rarefy is set to FALSE, compositional==TRUE, using relative abundances")
            x <- transformCounts(x, method = "relative", assay.type = assay.type)
        } else if (rarefy == TRUE && compositional == TRUE) {
            stop("Both rarefy and compositional cannot be TRUE simultaneously")
        }
        
        assay.data <- assay(x, assay.type)
        time <- colData(x)[, "Time.hr"]
        colnames(assay.data) <- time
        
        assay.data <- reshape2::melt(assay.data)
        
        ########################### Growth Metrics ############################
        grwt <- as_tibble(assay.data) %>%
            arrange(Var2) %>%
            group_by(Var1) %>%
            mutate(
                time_lag = Var2 - lag(Var2), 
                growth_diff = value - lag(value),
                growth_rate = value - lag(value) / lag(value),
                var_abund = value - lag(value) / time_lag
            )
        
        colnames(grwt)[colnames(grwt) == "Var1"] <- "OTU"
        colnames(grwt)[colnames(grwt) == "Var2"] <- "OTU"
        
        maxgrs <- grwt %>%
            summarize(max.growth = max(growth_diff, na.rm = T))
        grs.all <- merge(gwrt, maxgrs)
        grs.all <- grs.all %>%
            mutate(ismax = ifelse(growth_diff == max.growth, T, F))

        grs.all$OTU <- gsub("_", " ", grs.all$OTU)
        return(grs.all)
    }
)
        