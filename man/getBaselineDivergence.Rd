% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getBaselineDivergence.R
\name{getBaselineDivergence}
\alias{getBaselineDivergence}
\title{Beta diversity between the baseline and later time steps}
\usage{
getBaselineDivergence(
  x,
  group = NULL,
  time_field,
  name_divergence = "divergence_from_baseline",
  name_timedifference = "time_from_baseline",
  abund_values = "counts",
  FUN = vegan::vegdist,
  method = "bray",
  baseline_sample = NULL
)
}
\arguments{
\item{x}{A
\code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
object.}

\item{group}{a single character value for specifying which grouping
factor is used (name of a \code{colData} field). Optional.}

\item{time_field}{a single character value, specifying the name of the
time series field in \code{colData}.}

\item{name_divergence}{a column vector showing beta diversity between samples
over n time intervals (default: \code{name_divergence = "time_divergence"})}

\item{name_timedifference}{field name for adding the time difference between
samples used to calculate beta diversity
(default: \code{name_timedifference = "time_difference"})}

\item{abund_values}{character indicating which assay values are used in
the dissimilarity estimation (default: \code{abund_values = "counts"})}

\item{FUN}{a \code{function} for dissimilarity calculation. The function must
expect the input matrix as its first argument. With rows as samples
and columns as features. By default, \code{FUN} is
\code{vegan::vegdist}.}

\item{method}{a method that is used to calculate the distance. Method is
passed to the function that is specified by \code{FUN}. By default,
\code{method} is \code{"bray"}.}

\item{baseline_sample}{Optional. A character vector specifying the baseline sample(s) to be used. If the
"group" argument is given, this must be a named vector; one element per group.}
}
\value{
a
\code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
or
\code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
containing the sample dissimilarity and corresponding time difference between
samples (across n time steps), within each level of the grouping factor.
}
\description{
Calculates sample dissimilarity between the given baseline and other
time points, optionally within a group (subject, reaction chamber, or
similar). The corresponding time difference is returned as well.
The method operates on \code{SummarizedExperiment} objects, and the results
are stored in \code{colData}.
}
\details{
The group argument allows calculating divergence per group. Otherwise, this is done across all samples at once.

The baseline sample/s always need to belong to the data object i.e. they can be merged into it before
applying this function. The reason is that they need to have comparable sample data, at least some time point
information for calculating time differences w.r.t. baseline sample.

The baseline time point is by default defined as the smallest time point (per group). Alternatively,
the user can provide the baseline vector, or a list of baseline vectors per group (named list per group).
}
\examples{
#library(miaTime)
library(TreeSummarizedExperiment)
library(dplyr)

data(hitchip1006)
tse <- mia::transformSamples(hitchip1006, method = "relabundance")

# Subset to speed up example
tse <- tse[, colData(tse)$subject \%in\% c("900", "934", "843", "875")]

tse2 <- getBaselineDivergence(tse,
                              group = "subject",
                              time_field = "time",
                              name_divergence = "divergence_from_baseline",
                              name_timedifference = "time_from_baseline",
                              abund_values="relabundance",
                              FUN = vegan::vegdist,
                              method="bray")

tse2 <- getBaselineDivergence(tse,
                              baseline_sample = "Sample-875",
                              group = "subject",
                              time_field = "time",
                              name_divergence = "divergence_from_baseline",
                              name_timedifference = "time_from_baseline",
                              abund_values="relabundance",
                              FUN = vegan::vegdist,
                              method="bray")

}
