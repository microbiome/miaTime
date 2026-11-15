#' @name
#' getVariability
#'
#' @export
#'
#' @title
#' Variability of features
#'
#' @description
#' Calculates indices that describe how the abundance of each feature varies
#' over the samples, for example over time within each subject: how much
#' (variance, coefficient of variation), in what form (skewness, kurtosis,
#' bimodality) and how persistently (stability). These are the feature-level
#' counterparts of the sample-level indices of
#' \code{\link[mia:addAlpha]{getAlpha()}}.
#'
#' @details
#' The available indices are:
#' \itemize{
#'   \item \code{"variance"}: sample variance.
#'   \item \code{"cv"}: coefficient of variation, the standard deviation
#'   divided by the mean.
#'   \item \code{"skewness"}: moment coefficient of skewness.
#'   \item \code{"kurtosis"}: moment coefficient of kurtosis (3 for a normal
#'   distribution).
#'   \item \code{"bimodality"}: see
#'   \code{\link[miaTime:getBimodality]{getBimodality()}}.
#'   \item \code{"stability"}: see
#'   \code{\link[miaTime:getStability]{getStability()}}; requires
#'   \code{time.col}.
#' }
#' With \code{group}, each index is calculated within each group, for example
#' within each subject, and the columns are named \code{<index>_<group>}.
#'
#' @return
#' \code{getVariability} returns a \code{DataFrame} with one row per feature.
#' \code{addVariability} adds the results to \code{rowData(x)}.
#'
#' @param x \code{SummarizedExperiment}.
#'
#' @param assay.type \code{Character scalar}. Specifies the name of the assay
#' used in calculation. (Default: \code{"counts"})
#'
#' @param index \code{Character vector}. Specifies the indices to calculate.
#' (Default: \code{c("variance", "cv", "skewness", "kurtosis",
#' "bimodality")})
#'
#' @param group \code{NULL} or \code{Character scalar}. Specifies a column
#' from \code{colData(x)} that defines the groups, for example subjects,
#' within which the indices are calculated. (Default: \code{NULL})
#'
#' @param name \code{Character vector}. Specifies the column names for the
#' indices, one per index. (Default: \code{index})
#'
#' @param ... additional arguments passed to
#' \code{\link[miaTime:getBimodality]{getBimodality()}} and
#' \code{\link[miaTime:getStability]{getStability()}}, such as
#' \code{time.col}.
#'
#' @examples
#' library(miaTime)
#'
#' # A defined community of 16 species in three replicate bioreactors
#' data(minimalgut)
#' tse <- transformAssay(minimalgut, method = "relabundance")
#'
#' # Variability of each species over time, within each bioreactor
#' getVariability(tse, assay.type = "relabundance",
#'     index = c("cv", "stability"), group = "StudyIdentifier",
#'     time.col = "Time.hr")
#'
#' # Add the indices over all samples to rowData
#' tse <- addVariability(tse, assay.type = "relabundance")
#' rowData(tse)
#'
#' @seealso
#' \code{\link[mia:addAlpha]{getAlpha()}}
NULL

#' @rdname getVariability
#' @export
setMethod("getVariability", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = "counts", index = c(
            "variance", "cv", "skewness", "kurtosis", "bimodality"),
            group = NULL, ...){
        .get_variability(x, assay.type, index, name = index, group, ...)
    }
)

#' @rdname getVariability
#' @export
setMethod("addVariability", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = "counts", index = c(
            "variance", "cv", "skewness", "kurtosis", "bimodality"),
            name = index, group = NULL, ...){
        x <- .check_and_get_altExp(x, ...)
        res <- .get_variability(x, assay.type, index, name, group, ...)
        res <- as.list(res)
        x <- .add_values_to_colData(x, unname(res), names(res), MARGIN = 1L)
        return(x)
    }
)

################################ HELP FUNCTIONS ################################

# Indices that are calculated from the abundances of a feature alone.
.variability_functions <- list(
    variance = function(v) var(v),
    cv = function(v) sd(v) / mean(v),
    skewness = function(v) .calculate_skewness(v),
    kurtosis = function(v) .calculate_kurtosis(v)
)

.get_variability <- function(x, assay.type, index, name, group, ...){
    supported <- c(names(.variability_functions), "bimodality", "stability")
    temp <- .check_input(index, list("character vector"), supported)
    temp <- .check_input(name, list("character vector"),
        length = length(index))
    if( anyDuplicated(index) ){
        stop("'index' must not include duplicates.", call. = FALSE)
    }
    temp <- .check_input(
        group, list(NULL, "character scalar"), colnames(colData(x)))
    x <- .check_and_get_altExp(x, ...)
    .check_assay_present(assay.type, x)
    res <- mapply(function(i, nm){
        switch(i,
            bimodality = getBimodality(
                x, assay.type = assay.type, group = group, name = nm, ...),
            stability = getStability(
                x, assay.type = assay.type, group = group, name = nm, ...),
            .summarize_features(
                x, assay.type, .variability_functions[[i]], nm, group))
    }, index, name, SIMPLIFY = FALSE)
    res <- do.call(cbind, unname(res))
    return(res)
}

# Applies FUN to the abundances of each feature, within each group.
#' @importFrom dplyr group_by across all_of summarize
#' @importFrom tidyr pivot_wider
.summarize_features <- function(x, assay.type, FUN, name, group){
    .data <- value <- NULL
    df <- meltSE(x, assay.type = assay.type, add.col = group)
    res <- df |>
        group_by(across(all_of(c(group, "FeatureID")))) |>
        summarize(value = FUN(.data[[assay.type]]), .groups = "drop")
    if( !is.null(group) ){
        res <- pivot_wider(res, names_from = all_of(group),
            values_from = value)
    }
    res <- DataFrame(res, check.names = FALSE)
    rownames(res) <- res[["FeatureID"]]
    res[["FeatureID"]] <- NULL
    res <- res[match(rownames(x), rownames(res)), , drop = FALSE]
    if( ncol(res) == 1L ){
        colnames(res) <- name
    } else{
        colnames(res) <- paste0(name, "_", colnames(res))
    }
    return(res)
}
