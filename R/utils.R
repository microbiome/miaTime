################################################################################
# internal methods loaded from other packages

.check_altExp_present <- mia:::.check_altExp_present
.get_mat_from_sce <- scater:::.get_mat_from_sce

################################################################################
# internal wrappers for getter/setter

#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom S4Vectors DataFrame
.add_values_to_colData <- function(x, values, name){
    # converts each value:name pair into a DataFrame
    values <- mapply(
        function(value, n){
            value <- DataFrame(value)
            colnames(value)[1L] <- n
            if(ncol(value) > 1L){
                i <- seq.int(2,ncol(value))
                colnames(value)[i] <- paste0(n,"_",colnames(value)[i])
            }
            value
        },
        values,
        name)

    values <- do.call(cbind, values)

    # check for duplicated values
    f <- colnames(colData(x)) %in% colnames(values)
    if(any(f)) {
        warning("The following values are already present in `colData` and ",
                "will be overwritten: '",
                paste(colnames(colData(x))[f], collapse = "', '"),
                "'. Consider using the 'name' argument(s) to specify alternative ",
                "names.",
                call. = FALSE)
    }
    # keep only unique values
    colData(x) <- cbind(colData(x)[!f], values)

    x
}
