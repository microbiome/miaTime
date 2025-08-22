#' @name
#' getSurvival
#' 
#' @export
#' 
#' @title
#' Survival analysis
#' 
#' @description
#' Fit a (penalized) Cox proportional hazards model on microbiome data contained 
#' in a SummarizedExperiment object. 
#' Data transformations (e.g. pairwise log-ratios) should be handled upstream 
#' (e.g. with mia::transformAssay). This function focuses on the statistical 
#' model.
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
#' @param ... additional arguments.
#' \itemize{
#'   \item \code{penalized}: \code{Logical}. If TRUE, fit penalized Cox 
#'   regression using \code{glmnet}. If FALSE, fit standard Cox model using 
#'   \code{survival::coxph}. (Default: \code{TRUE})
#'
#'   \item \code{lambda}: \code{Character or numeric}. Penalization parameter 
#'   passed to \code{\link[glmnet]{cv.glmnet}}. Use \code{"lambda.1se"}, 
#'   \code{"lambda.min"}, or a numeric value. (Default: \code{"lambda.1se"})
#'
#'   \item \code{alpha}: \code{Numeric scalar}. Elastic net mixing parameter: 
#'   1 = Lasso, 0 = Ridge. (Default: \code{0.9})
#'
#'   \item \code{nfolds}: \code{Integer scalar}. Number of cross-validation 
#'   folds for \code{cv.glmnet}. (Default: \code{10})
#' 
#'   \item \code{nvar}: \code{Integer scalar}. Optional. Maximum number of 
#'   variables (log-ratios) to include in the model. (Default: \code{NULL})
#' 
#'   \item \code{coef_threshold}: \code{Numeric scalar}. Minimum absolute value 
#'   for a coefficient to be included in the final model. (Default: \code{0})
#' }  
#' 
#' @inheritParams addBaselineDivergence
#'
#' @return A list with model summaries:
#' \itemize{
#'   \item \code{coefficients}: estimated model coefficients
#'   \item \code{selected.features}: nonzero features in penalized model
#'   \item \code{risk.scores}: predicted risk scores
#'   \item \code{cindex.apparent}: apparent concordance index
#'   \item \code{cv.cindex.mean}: mean cross-validated C-index (if penalized)
#'   \item \code{cv.cindex.sd}: SD of cross-validated C-index (if penalized)
#'   \item \code{fit}: fitted model object (coxph or cv.glmnet)
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
#' @importFrom survival Surv coxph
#' @importFrom glmnet cv.glmnet glmnet
setMethod("getSurvival", signature(x = "SummarizedExperiment"),
    function(
        x,
        time.col,
        status.col,
        assay.type = "counts",
        covar.cols = NULL,
        ...
          ) {
        # input checks
        x <- .check_data_for_survival(x, time.col, 
            status.col, covar.cols, assay.type, ...)
        
        # extract data
        dat <- .get_data_for_survival(x, time.col, 
            status.col, covar.cols, assay.type)
        
        # fit survival model
        res <- .calc_survival(
            mat = dat$mat,
            time = dat$time,
            status = dat$status,
            covar = dat$covar,
            ...
        )
              
        return(res)
          }
)

############################# internal helpers #################################

.check_data_for_survival <- function(x, time.col, 
    status.col, covar.cols, assay.type, ...) 
{
    x <- .check_and_get_altExp(x, ...)
    .check_input(time.col, list("character scalar"), colnames(colData(x)))
    .check_input(status.col, list("character scalar"), colnames(colData(x)))
    .check_assay_present(assay.type, x)
    
    if (!is.numeric(x[[time.col]])) {
        stop("'time.col' must be numeric", call. = FALSE)
    }
    if (!is.logical(x[[status.col]]) && !is.numeric(x[[status.col]])) {
        stop("'status.col' must be numeric (0/1) or logical", call. = FALSE)
    }
    if (!is.null(covar.cols)) {
        .check_input(covar.cols, list("character vector"), colnames(colData(x)))
    }
    
    return(x)
}

.get_data_for_survival <- function(x, time.col, 
    status.col, covar.cols, assay.type) 
{
    mat <- t(assay(x, assay.type))
    time <- as.numeric(x[[time.col]])
    status <- as.numeric(x[[status.col]])
    
    covar <- NULL
    if (!is.null(covar.cols)) {
        covar <- as.data.frame(colData(x)[, covar.cols, drop = FALSE])
    }
    
    list(
        mat = mat,
        time = time,
        status = status,
        covar = covar
    )
}

.calc_survival <- function(mat, time, status, covar = NULL,
    penalized = TRUE, alpha = 0.9, nfolds = 10,
    lambda = "lambda.1se", nvar = NULL,
    coef_threshold = 0, ...) 
{
    y <- survival::Surv(time, status)
    
    if (penalized) {
        ## --- penalized Cox with glmnet ---
        if (is.null(covar)) {
            cvfit <- glmnet::cv.glmnet(
                x = mat,
                y = y,
                family = "cox",
                type.measure = "C",
                alpha = alpha,
                nfolds = nfolds,
                keep = TRUE
            )
            offset <- NULL
        } else {
            ## fit Cox model for covariates -> use offset
            df0 <- data.frame(time = time, status = status, covar)
            model0 <- survival::coxph(Surv(time, status) ~ ., data = df0)
            offset <- predict(model0)
            cvfit <- glmnet::cv.glmnet(
                x = mat,
                y = y,
                family = "cox",
                type.measure = "C",
                alpha = alpha,
                nfolds = nfolds,
                keep = TRUE,
                offset = offset
            )
        }
        
        ## --- nvar cutoff handling ---
        if (!is.null(nvar)) {
            rowlasso <- max(which(cvfit$glmnet.fit$df <= nvar))
            lambda <- cvfit$glmnet.fit$lambda[rowlasso]
        }
        
        ## --- lambda selection ---
        if (is.character(lambda)) {
            if (lambda == "lambda.1se") lambda.value <- cvfit$lambda.1se
            if (lambda == "lambda.min") lambda.value <- cvfit$lambda.min
        } else {
            lambda.value <- lambda
        }
        
        ## --- coefficients ---
        coefs <- as.vector(coef(cvfit, s = lambda.value))
        names(coefs) <- rownames(coef(cvfit, s = lambda.value))
        selected <- which(abs(coefs) > coef_threshold)
        coef.filtered <- coefs[selected]
        
        ## normalize coefficients
        if (length(coef.filtered) > 0) {
            coef.filtered <- 2 * coef.filtered / sum(abs(coef.filtered))
        }
        
        ## --- predictions ---
        if (is.null(offset)) {
            risk.scores <- as.numeric(predict(cvfit, mat, s = lambda.value))
        } else {
            risk.scores <- as.numeric(predict(cvfit, mat, s = lambda.value,
                                              newoffset = offset))
        }
        
        ## --- apparent C-index ---
        if (length(coef.filtered) == 0) {
            Cindex_signature <- 0.5
        } else {
            Cindex_signature <- glmnet::Cindex(pred = risk.scores, y)
        }
        
        ## --- CV stats at chosen lambda ---
        idrow <- max(which(cvfit$glmnet.fit$lambda >= lambda.value))
        mcvCindex <- cvfit$cvm[idrow]
        sdcvCindex <- cvfit$cvsd[idrow]
        
        out <- list(
            coefficients = coef.filtered,
            selected.features = names(coef.filtered),
            risk.scores = risk.scores,
            cindex.apparent = Cindex_signature,
            cv.cindex.mean = mcvCindex,
            cv.cindex.sd = sdcvCindex,
            fit = cvfit
        )
        
    } else {
        ## --- standard CoxPH ---
        df <- data.frame(time = time, status = status, mat, covar)
        formula <- as.formula(paste("Surv(time, status) ~ ."))
        fit <- survival::coxph(formula, data = df)
        
        ## raw coefficients
        coefs <- coef(fit)
        
        ## apply coef_threshold
        selected <- which(abs(coefs) > coef_threshold)
        coef.filtered <- coefs[selected]
        
        ## normalize
        if (length(coef.filtered) > 0) {
            coef.filtered <- 2 * coef.filtered / sum(abs(coef.filtered))
        }
        
        risk.scores <- predict(fit, type = "lp")
        
        ## apparent C-index
        if (length(coef.filtered) == 0) {
            Cindex_signature <- 0.5
        } else {
            Cindex_signature <- summary(fit)$concordance[1]
        }
        
        out <- list(
            coefficients = coef.filtered,
            selected.features = names(coef.filtered),
            risk.scores = as.numeric(risk.scores),
            cindex.apparent = Cindex_signature,
            cv.cindex.mean = NULL,
            cv.cindex.sd = NULL,
            fit = fit
        )
    }
    
    return(out)
}
