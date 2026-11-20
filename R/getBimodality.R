#' @name
#' getBimodality
#'
#' @title
#' Calculate coefficient of bimodality.
#'
#' @description
#' This function calculates coefficient of bimodality for each taxa.
#'
#' @details
#' This function calculates coefficient of bimodality for each taxa. If the
#' dataset includes grouping, for instance, individual systems or patients,
#' the coefficient is calculated for each group separately.
#'
#' The coefficient of bimodality measures whether a taxon has bimodal abundance.
#' For instance, certain taxon can be high-abundant in some timepoints while
#' in others it might be rare. The coefficient can help to determine these
#' taxa.
#'
#' The coefficient of bimodality (b) is defined as follows:
#'
#' \deqn{b = \frac{1+skewness^{2}}{kurtosis+3}}
#'
#' where skewness is calculated as follows
#'
#' \deqn{skewness = \frac{\sum_{i=1}^{n}(x_{i}-\overline{x})^{3}/n}{
#' \sum_{i=1}^{n}(x_{i}-\overline{x})^{2}/n]^{3/2}}}
#'
#' and kurtosis as follows
#'
#' \deqn{kurtosis = \frac{\sum_{i=1}^{n}(x_{i}-\overline{x})^{4}/n}{
#' (\sum_{i=1}^{n}(x_{i}-\overline{x})^{2}/n)^{2}}}
#'
#' The coefficient ranges from 0-1, where 1 means high bimodality. The
#' coefficient was introduced in the paper Shade A et al. 2014, where they used
#' bimodality and abundance to determine conditionally rare taxa (CRT).
#'
#' With \code{method = "potential"}, the score comes from potential analysis
#' (Livina et al. 2010), repeated on bootstrap resamples as in Lahti et al.
#' (2014). In each resample, the abundances of a feature are smoothed with a
#' Gaussian kernel density, and the maxima of the density that are high and
#' well separated by minima count as modes. The score is the share of
#' resamples with more than one mode, from 0 (unimodal) to 1. Missing values
#' are ignored. Small random noise separates tied values, so set a seed for
#' reproducible scores.
#'
#' @references
#'
#' Shade A, et al. (2014)
#' Conditionally Rare Taxa Disproportionately Contribute to Temporal Changes in
#' Microbial Diversity. doi: 10.1128/mbio.01371-14
#'
#' Livina V, et al. (2010)
#' Potential analysis reveals changing number of climate states during the
#' last 60 kyr. Climate of the Past 6, 77-82.
#'
#' Lahti L, et al. (2014)
#' Tipping elements in the human intestinal ecosystem.
#' Nature Communications 5, 4344.
#'
#' @return
#' \code{DataFrame} or \code{x} with results added to its \code{rowData}.
#'
#' @inheritParams addBaselineDivergence
#'
#' @param method \code{Character scalar}. The score: the coefficient of
#' bimodality (\code{"coefficient"}) or the multimodality score of potential
#' analysis (\code{"potential"}). (Default: \code{"coefficient"})
#'
#' @param name \code{Character scalar}. Specifies a column name for storing
#' bimodality results. (Default: \code{"bimodality"})
#'
#' @param ... additional arguments.
#' \itemize{
#'   \item \code{group}: \code{Character scalar}. Specifies a name of the column
#'   from \code{colData} that identifies the grouping of the samples.
#'   (Default: \code{NULL})
#'   \item \code{peak.threshold}: \code{Numeric scalar}. For
#'   \code{method = "potential"}, how much higher than its neighbouring
#'   minima a maximum must be to count as a mode, in units of the kernel
#'   height. (Default: \code{1})
#'   \item \code{bw.adjust}: \code{Numeric scalar}. For
#'   \code{method = "potential"}, multiplies the bandwidth of the density.
#'   (Default: \code{1})
#'   \item \code{bs.iter}: \code{Integer scalar}. For
#'   \code{method = "potential"}, the number of bootstrap resamples.
#'   (Default: \code{100})
#'   \item \code{min.density}: \code{Numeric scalar}. For
#'   \code{method = "potential"}, the lowest density of a mode, in units of
#'   the kernel height. (Default: \code{1})
#' }
#'
#' @examples
#' library(miaTime)
#'
#' data(SilvermanAGutData)
#' tse <- SilvermanAGutData
#'
#' # In this example, we are only interested on vessel 1, in samples with reads.
#' tse <- tse[, tse[["Vessel"]] == 1 & colSums(assay(tse, "counts")) > 0]
#' tse <- transformAssay(tse, method = "relabundance")
#'
#' # Calculate bimodality
#' b <- getBimodality(tse)
#' # Determine taxa with high bimodality
#' bimodal_taxa <- names(b)[ which(b[[1]] > 0.95) ]
#'
#' # Determine taxa with abundance > 0.5%
#' abundant <- getAbundant(
#'     tse, assay.type = "relabundance", abundant.th = 0.5/100)
#'
#' # The detected CRT
#' crt <- intersect(bimodal_taxa, abundant)
#' head(crt)
#'
#' # Multimodality from potential analysis, for the abundant taxa
#' set.seed(1)
#' p <- getBimodality(tse[abundant, ], assay.type = "relabundance",
#'     method = "potential", bs.iter = 20)
#' head(p[order(-p[[1]]), , drop = FALSE])
#'
#' @seealso
#' \code{\link[mia:getAbundant]{mia::getConditionallyLowAbundant()}}
#'
NULL

#' @rdname getBimodality
#' @export
setMethod("addBimodality", signature = c(x = "SummarizedExperiment"),
    function(x, name = "bimodality", ...){
        .check_input(name, "character scalar")
        x <- .check_and_get_altExp(x, ...)
        args <- c(
            list(x = x, name = name),
            list(...)[!names(list(...)) %in% c("altexp")])
        res <- do.call(getBimodality, args) |> as.list()
        x <- .add_values_to_colData(x, unname(res), names(res), MARGIN = 1L)
        return(x)
    }
)

#' @rdname getBimodality
#' @export
setMethod("getBimodality", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = "counts",
            method = c("coefficient", "potential"), ...){
        x <- .check_and_get_altExp(x, ...)
        .check_assay_present(assay.type, x)
        method <- match.arg(method)
        res <- .calculate_bimodality(x, assay.type, method, ...)
        return(res)
    }
)

################################ HELP FUNCTIONS ################################

# This function calculates bimodality for each group for eacg feature.
#' @importFrom dplyr group_by across all_of summarize
#' @importFrom tidyr pivot_wider
.calculate_bimodality <- function(
        tse, assay.type, method = "coefficient", name = "bimodality",
        group = NULL, ...){
    # This following line is to suppress "no visible binding for" messages
    # in cmdcheck
    .data <- b <- NULL

    .check_input(name, "character scalar")
    .check_input(group, c("character scalar", "NULL"), colnames(colData(tse)))
    # Get data as long format
    df <- meltSE(tse, assay.type = assay.type, add.col = group)
    # Calculate bimodality for each group and feature
    res <- df |>
        group_by(across(all_of(c(group, "FeatureID")))) |>
        summarize(
            b = switch(method,
                coefficient = .calculate_bimodality_coefficient(
                    .data[[assay.type]], ...),
                potential = .calculate_multimodality(.data[[assay.type]], ...)),
            .groups = "drop"
            )
    # If grouping was used, put data to wide format so that columns include
    # groups and rows features.
    if( !is.null(group) ){
        res <- res |>
            pivot_wider(
                names_from = group,
                values_from = b
            )
    }
    # Put into DataFrame format and wrangle rownames
    res <- DataFrame(res, check.names = FALSE)
    rownames(res) <- res[["FeatureID"]]
    res[["FeatureID"]] <- NULL
    # Make sure that the order of data reflects the original data.
    res <- res[match(rownames(tse), rownames(res)), , drop = FALSE]
    # Rename columns by user specified name
    if( ncol(res) == 1L ){
        colnames(res) <- name
    } else{
        colnames(res) <- paste0(name, "_", colnames(res))
    }
    return(res)
}

# Function to calculate bimodality coefficient for single vector/taxon
.calculate_bimodality_coefficient <- function(x, ...) {
    skewness <- .calculate_skewness(x, ...)
    kurtosis <- .calculate_kurtosis(x, ...)
    b <- (1 + skewness^2) / (kurtosis + 3)
    return(b)
}

# Function to calculate skewness for single vector/taxon
.calculate_skewness <- function(x, ...) {
    mean_x <- mean(x, ...)
    n <- length(x)
    numerator <- sum((x - mean_x)^3) / n
    denominator <- (sum((x - mean_x)^2) / n)^(3/2)
    skewness <- numerator / denominator
    return(skewness)
}

# Function to calculate kurtosis for single vector/taxon
.calculate_kurtosis <- function(x, ...) {
    mean_x <- mean(x, ...)
    n <- length(x)
    numerator <- sum((x - mean_x)^4) / n
    denominator <- (sum((x - mean_x)^2) / n)^2
    kurtosis <- numerator / denominator
    return(kurtosis)
}

# Multimodality score of one feature: the share of bootstrap resamples with
# more than one mode in the potential analysis (as in microbiome).
#' @importFrom stats rnorm
.calculate_multimodality <- function(
        x, peak.threshold = 1, bw.adjust = 1, bs.iter = 100,
        min.density = 1, ...){
    x <- x[!is.na(x)]
    if( length(unique(x)) <= 1L ){
        return(0)
    }
    # Small noise separates tied values
    x <- x + rnorm(length(x), sd = sd(x) / 100)
    n_modes <- vapply(seq_len(bs.iter), function(i){
        xb <- x[sample(length(x), replace = TRUE)]
        length(.potential_maxima(xb, peak.threshold, bw.adjust, min.density))
    }, numeric(1))
    return(1 - mean(n_modes <= 1))
}

# Modes of the kernel density of x: maxima that are dense enough and well
# separated by minima.
#' @importFrom stats density dnorm
.potential_maxima <- function(x, peak.threshold, bw.adjust, min.density){
    de <- density(x, bw = "nrd", adjust = bw.adjust, kernel = "gaussian",
        n = floor(0.2 * length(x)), from = min(x), to = max(x), cut = 3)
    f <- de[["y"]] / sum(diff(de[["x"]][seq_len(2)]) * de[["y"]])
    height <- dnorm(0, sd = de[["bw"]]) / length(f)
    maxima <- .find_maxima(f)
    maxima <- maxima[f[maxima] >= min.density * height]
    minima <- .between_minima(f, maxima, .find_maxima(-f))
    # Drop maxima that rise less than the threshold above a neighbouring
    # minimum, together with the higher of their two minima
    drop_max <- logical(length(maxima))
    drop_min <- logical(length(minima))
    for( j in seq_along(maxima) ){
        right <- minima[minima > maxima[j] & !drop_min]
        left <- minima[minima < maxima[j] & !drop_min]
        i1 <- if( length(right) ) min(right) else if( length(left) )
            max(left)
        i2 <- if( length(left) ) max(left) else i1
        if( !is.null(i1) &&
                min(f[maxima[j]] - f[i1], f[maxima[j]] - f[i2]) <
                peak.threshold * height ){
            drop_max[j] <- TRUE
            drop_min[minima == if( f[i1] > f[i2] ) i1 else i2] <- TRUE
        }
    }
    return(maxima[!drop_max])
}

# Indices of local maxima of f; a flat top gives its middle index.
.find_maxima <- function(f){
    d <- sign(diff(c(-Inf, f, -Inf)))
    # Plateaus take the direction of the next change
    for( i in rev(which(d == 0)) ){
        d[i] <- d[i + 1L]
    }
    ends <- which(diff(d) < 0)
    starts <- vapply(ends, function(e){
        s <- e
        while( s > 1L && f[s - 1L] == f[e] ){
            s <- s - 1L
        }
        s
    }, numeric(1))
    return(round((starts + ends) / 2))
}

# The lowest minimum between each pair of consecutive maxima.
.between_minima <- function(f, maxima, minima){
    if( length(maxima) < 2L ){
        return(integer(0))
    }
    res <- lapply(seq_along(maxima)[-1L], function(i){
        m <- minima[minima > maxima[i - 1L] & minima < maxima[i]]
        if( length(m) ) round(mean(m[f[m] == min(f[m])]))
    })
    return(unlist(res))
}
