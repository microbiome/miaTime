#' @name
#' getSurvival
#' 
#' @export
#' 
#' @title
#' Survival analysis
#' 
#' @description
#' Survival modeling wrapper for TreeSummarizedExperiment using coda_coxnet
#' 
#' @details
#' This function extracts compositional count data and survival metadata from a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object and applies penalized Cox regression via \code{\link[coda4microbiome]{coda_coxnet}}, 
#' using pairwise log-ratio transformations.
#'
#'
#' @param time.col \code{Character scalar}. Column name in \code{colData(x)} 
#' representing time to event or follow-up time. Must be numeric.
#'
#' @param status.col \code{Character scalar}. Column name in \code{colData(x)} 
#' representing event occurrence. Accepts numeric (0/1) or logical (FALSE/TRUE) 
#' values.
#'
#' @param covar.cols \code{Character vector}. Optional. Specifies covariate 
#' columns in \code{colData(x)} to adjust for in the survival model. 
#' (Default: \code{NULL})
#'
#' @param lambda \code{Character or numeric}. Penalization parameter passed to 
#' \code{\link[glmnet]{cv.glmnet}}. Use \code{"lambda.1se"}, 
#' \code{"lambda.min"}, or a numeric value. (Default: \code{"lambda.1se"})
#'
#' @param nvar \code{Integer scalar}. Optional. Maximum number of variables 
#' (log-ratios) to include in the model. (Default: \code{NULL})
#'
#' @param alpha \code{Numeric scalar}. Elastic net mixing parameter: 1 = Lasso, 
#' 0 = Ridge. (Default: \code{0.9})
#'
#' @param nfolds \code{Integer scalar}. Number of cross-validation folds for 
#' \code{cv.glmnet}. (Default: \code{10})
#'
#' @param showPlots \code{Logical}. If \code{TRUE}, generates plots for 
#' cross-validation, risk scores, and selected signature. (Default: \code{TRUE})
#'
#' @param coef_threshold \code{Numeric scalar}. Minimum absolute value for a 
#' coefficient to be included in the final model. (Default: \code{0})
#' 
#' @inheritParams addBaselineDivergence
#'
#' @return A named \code{list} containing:
#' \itemize{
#'   \item \code{taxa.num}: indices of selected taxa
#'   \item \code{taxa.name}: names of selected taxa
#'   \item \code{log-contrast coefficients}: coefficients of selected taxa in 
#'   the final log-contrast model
#'   \item \code{risk.score}: predicted risk scores
#'   \item \code{apparent Cindex}: concordance index on training data
#'   \item \code{mean cv-Cindex}: mean cross-validated concordance index
#'   \item \code{sd cv-Cindex}: standard deviation of cross-validated 
#'   concordance index
#'   \item \code{risk score plot}: ggplot object showing risk score vs survival
#'   \item \code{signature plot}: ggplot object of selected taxa and coefficients
#' }
#'
#' @seealso \code{\link[coda4microbiome]{coda_coxnet}}, 
#' \code{\link[glmnet]{cv.glmnet}}
#' 
#' @references
#' Meritxell Pujolassos , Antoni Susín , M.Luz Calle (2024). 
#' \emph{Microbiome compositional data analysis for survival studies }
#' NAR Genomics and Bioinformatics, 6(2), lqae038. 
#' \doi{10.1093/nargab/lqae038}
#'
#' @examples
#' # data(SurvivalData)
#' # tse <- SurvivalData
#' # getSurvival(tse, time.col = "T1Dweek", status.col = "T1D", 
#' # covar.cols = c("Sex", "Antibiotics"))
#'
NULL

#' @rdname getSurvival
#' @export
#' @importFrom coda4microbiome coda_coxnet
setMethod("getSurvival", signature(x = "SummarizedExperiment"),
    function(
        x,
        time.col,
        status.col,
        assay.type = "counts",
        covar.cols = NULL,
        lambda = "lambda.1se",
        nvar = NULL,
        alpha = 0.9,
        nfolds = 10,
        showPlots = TRUE,
        coef_threshold = 0,
        ...
          ) {
              ############################# INPUT CHECK ########################
              x <- .check_and_get_altExp(x, ...)
              
              # Check required colData fields
              .check_input(time.col, list("character scalar"), 
                           colnames(colData(x)))
              .check_input(status.col, list("character scalar"), 
                           colnames(colData(x)))
              
              if (!is.numeric(x[[time.col]])) {
                  stop("'time.col' must specify a numeric column from colData(x)", 
                       call. = FALSE)
              }
              if (!is.logical(x[[status.col]]) && !is.numeric(x[[status.col]])) {
                  stop("'status.col' must be numeric (0/1) or logical (FALSE/TRUE)", 
                       call. = FALSE)
              }
              
              .check_assay_present(assay.type, x)
              
              if (!is.null(covar.cols)) {
                  .check_input(covar.cols, list("character vector"), 
                               colnames(colData(x)))
              }
              
              ########################### INPUT CHECK END ######################
              
              # Extract matrix
              mat <- t(assay(x, assay.type))
              
              # Extract survival time and status
              time <- as.numeric(x[[time.col]])
              status <- as.numeric(x[[status.col]])
              
              # Extract covariates (if specified)
              covar <- NULL
              if (!is.null(covar.cols)) {
                  covar <- as.data.frame(colData(x)[, covar.cols, drop = FALSE])
              }
              
              # Run core model
              result <- coda_coxnet(
                  x = mat,
                  time = time,
                  status = status,
                  covar = covar,
                  lambda = lambda,
                  nvar = nvar,
                  alpha = alpha,
                  nfolds = nfolds,
                  showPlots = showPlots,
                  coef_threshold = coef_threshold
              )
              
              return(result)
          }
)