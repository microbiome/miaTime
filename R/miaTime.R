#' \code{miaTime} Package.
#'
#' \code{miaTime} implements time series methods in \link[mia]{mia} ecosystem.
#' @name miaTime-package
#' @aliases miaTime
#' @seealso
#' \link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment} class
"_PACKAGE"
NULL

#' @import mia
#' @import methods
#' @import S4Vectors
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import TreeSummarizedExperiment
NULL

#' @title HITChip Atlas with 1006 Western Adults
#' @description This data set contains genus-level microbiota profiling with
#' HITChip for 1006 western adults with no reported health complications,
#' reported in Lahti et al. (2014).
#' @name hitchip1006
#' @details The data is also available for download from the Data Dryad
#' \url{http://doi.org/10.5061/dryad.pk75d} and was initially released
#' as a phyloseq object under the name atlas1006 in \pkg{microbiome} R package.
#' The \pkg{miaTime} package provides this as
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object.
#' Some of the subjects also have short time series.
#' @docType data
#' @author Leo Lahti
#' @return Loads the data set in R.
#' @references
#' Lahti et al. Tipping elements of the human intestinal ecosystem.
#' Nature Communications 5:4344, 2014. \url{https://doi.org/10.1038/ncomms5344}.
#' @usage data(hitchip1006)
#' @format The data set in
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
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
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object.
#' @docType data
#' @author Doris Vandeputte
#' @return Loads the data set in R.
#' @references
#' Vandeputte, D., De Commer, L., Tito, R.Y. et al. Temporal variability
#' in quantitative human gut microbiome profiles and implications for clinical
#' research.
#' Nat Commun 12, 6740 (2021). https://doi.org/10.1038/s41467-021-27098-7
#' @usage data(temporalMicrobiome20)
#' @format The data set in
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' format.
#' @keywords data
#'
NULL

#' @title Human Gut Minimal Microbiome Profiling Data
#' @description This time-series data set contains absolute temporal abundances
#' for 16 human gut species that were assembled in an in vitro gut system. These
#' were subjected to a variety of disturbances over a period of 20 days. The
#' sample data includes measurements for Acetate, Butyrate, Propionate, and
#' optical density.
#' @name minimalgut
#' @details The data is available also on
#' \url{https://github.com/microsud/syncomR}
#' The data contains 183 samples from 3 in vitro gut system, 61 per bioreactor
#' collected three times daily.
#' The current implementation in \pkg{miaTime} provides this as
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object.
#' @docType data
#' @author Sudarshan A. Shetty
#' @return Loads the data set in R.
#' @references
#' Shetty, S.A., Kostopoulos, I., Geerlings, S.Y. et al. Dynamic metabolic
#' interactions and trophic roles of human gut microbes identified using a
#' minimal microbiome exhibiting ecological properties.
#' ISME, 16, 2144–2159 (2022). https://doi.org/10.1038/s41396-022-01255-2.
#' @usage data(minimalgut)
#' @format The data set in
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' format.
#' @keywords data
#'
NULL

#' @title SilvermanAGutData
#' @description
#' The SilvermanAGutData dataset contains 16S rRNA gene based high-throughput
#' profiling of 4 in vitro artificial gut models. The sampling was done hourly
#' and daily to capture sub-daily dynamics of microbial community originating
#' from human feces. The data consists of 413 taxa from 639 samples.
#' The data set can be used to investigate longitudinal dynamics of microbial
#' community in a controlled environment.
#'
#' Column metadata includes the days of sampling, vessel identifier, sampling
#' frequency pre-post challenge with Bacteroides ovatus.
#'
#' The row metadata of the microbiome data contains taxonomic information on
#' the Kingdom, Phylum, Class, Order, Family and Genus and Species level.
#'
#' The row tree consists of a phylogenetic tree build using sequence
#' information of 413 taxa.
#'
#' As reference sequences the ASV are provided.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use
#'
#' @name SilvermanAGutData
#' @docType data
#' @author Sudarshan A. Shetty and Felix G.M. Ernst
#' @return Loads the data set in R.
#' @references
#' Silveman J.D et al. (2018): Dynamic linear models guide design and analysis
#' of microbiota studies within artificial human guts.
#' Microbiome 6:202 \url{https://doi.org/10.1186/s40168-018-0584-3}
#' @usage data(SilvermanAGutData)
#' @format The data set in
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' format.
#' @keywords data
#'
NULL

#' @title Kumaraswamy2024
#' @description
#' The Kumaraswamy2024 includes microbiota and metabolite profiling data 
#' from 78 Indian individuals (40 males, 38 females).
#' 
#' The Indian subjects were grouped into four diet groups (~20 subjects per group), 
#' and fecal samples were collected across three seasonal time points.
#'
#' The microbiota profiling was performed using HITChip microarray analysis 
#' (in duplicate), qPCR (in triplicate with eight-point standard curves), and 
#' LC-HRMS and HPLC metabolite profiling with internal standards.
#'
#' Column metadata includes diet group assignment, sampling season, sex, BMI, 
#' age, and questionnaire-based lifestyle metadata.
#'
#' Quality control metrics include Pearson correlation (>0.98) for HITChip, 
#' qPCR assay efficiency (>0.99), and technical replicates for HPLC and qPCR.
#'
#' Data sources:
#' - Microbiota HITChip microarray data
#' - qPCR absolute abundance data
#' - Chemical profiling data (HPLC, LC-HRMS)
#' - Sample metadata (diet, lifestyle)
#'
#' Processed and raw data are available via:
#' - Zenodo (DOI: https://doi.org/10.5281/zenodo.14424024)
#' - NCBI-SRA (fermented foods 16S rRNA sequencing, accession: PRJNA1191989)
#'
#' @name Kumaraswamy2024
#' @docType data
#' @author Jeyaram, K., Lahti, L., Tims, S. et al
#' @return Loads the data set in R.
#' @references
#' Jeyaram, K., Lahti, L., Tims, S. et al. Fermented foods affect the seasonal 
#' stability of gut bacteria in an Indian rural population.
#' Nat Commun 16, 771 \url{https://doi.org/10.1038/s41467-025-56014-6}
#' @usage data(Kumaraswamy2024)
#' @format The data set in
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' format.
#' @keywords data
#'
NULL

#' @name crohn_survival
#' @title
#' Survival microbiome data from 150 individuals with longitudinal measurements
#'
#' @description
#' Simulated dataset based on a Crohn's disease microbiome study. The dataset
#' is right-censored and includes time-to-event data, representing either the
#' occurrence of disease or censoring. It contains 150 individuals and 48 taxa.
#' The dataset originates from the \pkg{coda4microbiome} package under the
#' name "data_survival".
#'
#' @details
#' Sample metadata includes the following information:
#'  - diagnosis: Indicates whether an individual has Crohn’s disease (CD) or
#'    is a control.
#'  - event: Binary variable indicating disease status (1 = Crohn’s Disease,
#'    0 = Control).
#'  - event_time: The time of event occurrence or censoring (unitless time).
#'
#' Taxa data do not include additional taxonomy information.
#'
#' @references
#' Calle ML et al. (2023) coda4microbiome: compositional data analysis for
#' microbiome cross-sectional and longitudinal studies. BMC Bioinformatics,
#' 24(82). \url{https://doi.org/10.1186/s12859-023-05205-3}
#'
#' Gevers D et al. (2018) The Treatment-Naive Microbiome in New-Onset Crohn’s
#' Disease. Cell Host & Microbe, 15(4).
#' #' \url{https://doi.org/10.1016/j.chom.2014.02.005}
#'
#' @docType data
#'
#' @return
#' Loads the dataset in R.
#'
#' @usage
#' data(crohn_survival)
#'
#' @format
#' The dataset is in the
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' format.
#'
#' @keywords
#' data
NULL
