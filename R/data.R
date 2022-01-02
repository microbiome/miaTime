#' @title HITChip Atlas with 1006 Western Adults
#' @description This data set contains genus-level microbiota profiling with
#' HITChip for 1006 western adults with no reported health complications,
#' reported in Lahti et al. (2014).
#' @name hitchip1006
#' @details The data is also available for download from the Data Dryad
#' \url{http://doi.org/10.5061/dryad.pk75d} and was initially released
#' as a phyloseq object under the name atlas1006 in \pkg{microbiome} R package.
#' The current implementation in \pkg{miaTime} provides this as
#' \code{\link{TreeSummarizedExperiment-class}} object.
#' Some of the subjects also have short time series.
#' @docType data
#' @author Leo Lahti
#' @return Loads the data set in R.
#' @references
#' Lahti et al. Tipping elements of the human intestinal ecosystem.
#' Nature Communications 5:4344, 2014. \url{https://doi.org/10.1038/ncomms5344}.
#' @usage data(hitchip1006)
#' @format The data set in \code{\link{TreeSummarizedExperiment-class}}
#' format.
#' @keywords data
#'
NULL

#' @title Gut Microbiome Profiling 20 Belgian Women
#' @description This time-series data set contains absolute temporal variation
#' for most major gut genera, as well as diversity indicators for both
#' relative and quantitative abundance profiles of healthy women living in an
#' industrialized country.
#' @name temporalMicrobiome20
#' @details The data is available also on
#' \url{https://www.nature.com/articles/s41467-021-27098-7#Sec43}
#' The data contains 713 fecal samples from 20 Belgian women collected
#' over six weeks.
#' The current implementation in \pkg{miaTime} provides this as
#' \code{\link{TreeSummarizedExperiment-class}} object.
#' @docType data
#' @author Doris Vandeputte
#' @return Loads the data set in R.
#' @references
#' Vandeputte, D., De Commer, L., Tito, R.Y. et al. Temporal variability
#' in quantitative human gut microbiome profiles and implications for clinical
#' research.
#' Nat Commun 12, 6740 (2021). https://doi.org/10.1038/s41467-021-27098-7
#' @usage data(temporalMicrobiome20)
#' @format The data set in \code{\link{TreeSummarizedExperiment-class}}
#' format.
#' @keywords data
#'
NULL
