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
#' in a SummarizedExperiment object. Data transformations (e.g. pairwise
#' log-ratios) should be handled upstream
#' (e.g. with \code{mia::transformAssay()}).
#'
#' @param time.col \code{Character scalar}. Column name in \code{colData(x)}
#' representing time to event or follow-up time. Must be numeric.
#'
#' @param status.col \code{Character scalar}. Column name in \code{colData(x)}
#' representing event occurrence. Accepts numeric (\code{0}/\code{1}) or
#' logical (\code{TRUE}/\code{FALSE}) values.
#'
#' @param col.var \code{Character vector}. Optional. Specifies covariate
#' columns in \code{colData(x)} to adjust for in the survival model.
#' (Default: \code{NULL})
#'
#' @param ... additional arguments.
#' \itemize{
#'   \item \code{penalized}: \code{Logical}. If \code{TRUE}, fit penalized Cox
#'   regression using \code{glmnet}. If \code{FALSE}, fit standard Cox model
#'   using \code{survival::coxph}. (Default: \code{TRUE})
#'
#'   \item \code{lambda}: \code{Character or numeric}. Penalization parameter
#'   passed to \code{\link[glmnet]{cv.glmnet}}. Use \code{"lambda.1se"},
#'   \code{"lambda.min"}, or a numeric value. (Default: \code{"lambda.1se"})
#'
#'   \item \code{alpha}: \code{Numeric scalar}. Elastic net mixing parameter
#'   that controls the balance between Lasso and Ridge regression:
#'   \code{alpha = 1} corresponds to Lasso,
#'   \code{alpha = 0} corresponds to Ridge.
#'   Values between 0 and 1 specify a combination of the two.
#'   (Default: \code{0.9})
#'
#'   \item \code{nfolds}: \code{Integer scalar}. Number of cross-validation
#'   folds for \code{cv.glmnet}. (Default: \code{10})
#'
#'   \item \code{nvar}: \code{Integer scalar}. Optional. Maximum number of
#'   variables (log-ratios) to include in the model. (Default: \code{NULL})
#'
#'   \item \code{coef.threshold}: \code{Numeric scalar}. Minimum absolute value
#'   for a coefficient to be included in the final model. (Default: \code{0})
#' }
#'
#' @inheritParams addBaselineDivergence
#'
#' @return A list with model summaries:
#' \itemize{
#'   \item \code{coefficients}: estimated model coefficients
#'   \item \code{risk_scores}: predicted risk scores
#'   \item \code{c_index}: apparent concordance index
#'   \item \code{c_index_cv_mean}: mean cross-validated C-index (if penalized)
#'   \item \code{c_index_cv_sd}: SD of cross-validated C-index (if penalized)
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
#' data(crohn_survival)
#' tse <- crohn_survival
#' tse <- transformAssay(tse, method = "relabundance")
#' fit <- getSurvival(
#'     tse, assay.type = "relabundance",
#'     time.col = "event_time", status.col = "event"
#' )
#'
NULL

#' @rdname getSurvival
#' @export
setMethod("addSurvival", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, status.col, name = "survival", ...){
        .check_input(name, "character scalar")
        x <- .check_and_get_altExp(x, ...)
        # Run analysis
        args <- c(
            list(x = x, time.col = time.col, status.col = status.col,
                name = name),
            list(...)[!names(list(...)) %in% c("altexp")])
        res <- do.call(getSurvival, args)
        # Add results to metadata
        x <- .add_values_to_metadata(x, name, res, ...)
        return(x)
    }
)

#' @rdname getSurvival
#' @export
setMethod("getSurvival", signature(x = "SummarizedExperiment"),
    function(x, time.col, status.col, assay.type = "counts", col.var = NULL,
        ...){
        # Input checks
        x <- .check_data_for_survival(
            x, time.col, status.col, col.var, assay.type, ...)
        # Extract data
        args <- .get_data_for_survival(
            x, time.col, status.col, col.var, assay.type)
        args <- c(args, list(...))
        # Fit survival model
        res <- do.call(.calc_survival, args)
        return(res)
    }
)

############################# Internal helpers #################################

# Check input validity for survival analysis
#
# Ensures that the required columns for survival analysis are present in
# colData, are of the correct type, and that the requested assay exists.
.check_data_for_survival <- function(
        x, time.col, status.col, col.var, assay.type, ...){
    # Ensure we are working with the correct alternative experiment
    x <- .check_and_get_altExp(x, ...)
    # Check that 'time.col' exists in colData and is a character scalar
    .check_input(time.col, list("character scalar"), colnames(colData(x)))
    # Check that 'status.col' exists in colData and is a character scalar
    .check_input(status.col, list("character scalar"), colnames(colData(x)))
    # Check that the requested assay is present in the object
    .check_assay_present(assay.type, x)
    # Verify that the time column contains numeric values
    if( !is.numeric(x[[time.col]]) ){
        stop("'time.col' must be numeric.", call. = FALSE)
    }
    # Verify that the status column is either logical or numeric (0/1)
    if( !is.logical(x[[status.col]]) && !is.numeric(x[[status.col]]) ){
        stop("'status.col' must be numeric (0/1) or logical.", call. = FALSE)
    }
    # If a grouping variable is provided, check it exists in colData
    if( !is.null(col.var) ){
        .check_input(col.var, list("character vector"), colnames(colData(x)))
    }
    return(x)
}

# Extract and prepare data for survival analysis
#
# Retrieves the assay matrix, survival time, survival status,
# and optional covariates from the input object,
# formatted for downstream survival models.
.get_data_for_survival <- function(
        x, time.col, status.col, col.var, assay.type){
    # Extract assay data (samples as rows, features as columns)
    mat <- assay(x, assay.type) |> t()
    # Extract survival time column and ensure numeric
    time <- x[[time.col]] |> as.numeric()
    # Extract survival status column and ensure numeric (0/1)
    status <- x[[status.col]] |> as.numeric()
    # Extract optional covariates from colData
    if( !is.null(col.var) ){
        col.var <- colData(x)[, col.var, drop = FALSE] |> as.data.frame()
    }
    # Return prepared components in a list
    res <- list(
        mat = mat,
        time = time,
        status = status,
        col.var = col.var
    )
    return(res)
}


# Dispatcher to fit a survival model
#
# Chooses between a penalized Cox model (via glmnet) and a standard Cox
# proportional hazards model, based on the `penalized` argument.
.calc_survival <- function(mat, time, status, col.var, penalized = TRUE, ...){
    if( !.is_a_bool(penalized) ){
        stop("'penalized' must be TRUE or FALSE.", call. = FALSE)
    }
    FUN <- if( penalized ) .fit_penalized_cox else .fit_standard_cox
    res <- FUN(mat = mat, time = time, status = status, col.var = col.var, ...)
    return(res)
}

# Fit a penalized Cox proportional hazards model using glmnet
#
# Performs cross-validated elastic net penalized Cox regression. Optional
# covariates
# can be included via an offset.
#' @importFrom survival Surv coxph
#' @importFrom glmnet cv.glmnet glmnet
.fit_penalized_cox <- function(
        mat, time, status, col.var, alpha = 0.9, nfolds = 10,
        lambda = "lambda.1se", ...){
    # Create survival response object as glmnet requires a Surv object for Cox
    # regression
    y <- Surv(time, status)

    # Fit Cox model for covariates (if provided) and use as offset. This is done
    # to estimate how microbiome profile improves the prediction compared to
    # conventional clinical variables.
    offset <- NULL
    if( !is.null(col.var) ){
        # Fit standard Cox on covariates only
        df_covar <- data.frame(time = time, status = status, col.var)
        model_covar <- coxph(y ~ ., data = df_covar)
        # Compute the linear predictor (log-hazard ratio) from the
        # covariate-only Cox model.
        offset <- predict(model_covar, type = "lp")
    }

    # Fit penalized Cox model with cross-validation. Cross-validation identifies
    # the optimal penalty (lambda) while controlling overfitting. This model
    # uses microbial profile as predictor. If covariates where specified, the
    # model is built on top of the covariate-only model, i.e., to assess how
    # much microbes improve the prediction.
    fit <- cv.glmnet(
        x = mat, y = y, family = "cox", type.measure = "C",
        alpha = alpha, nfolds = nfolds, keep = TRUE,
        offset = offset
    )

    # Select lambda value based on user input and optional nvar constraint to
    # allow selecting lambda that balances model sparsity and predictive
    # performance
    lambda_value <- .lambda_selector(fit, lambda, ...)
    # Extract coefficients at chosen lambda
    coefs <- coef(fit, s = lambda_value)
    coefs <- setNames(as.vector(coefs), rownames(coefs))
    coefs <- .filter_and_normalize_coefs(coefs, ...)

    # Compute risk scores for all samples
    risk_scores <- predict(fit, mat, s = lambda_value, newoffset = offset) |>
        as.numeric()
    # Compute concordance index
    c_index <- .compute_c_index(risk_scores, y, penalized = TRUE)

    # Identify row in CV results corresponding to chosen lambda
    id_row <- which(fit[["glmnet.fit"]][["lambda"]] >= lambda_value) |> max()

    # Return all results as a list
    res <- list(
        coefficients = coefs,
        risk_scores = risk_scores,
        c_index = c_index,
        c_index_cv_mean = fit[["cvm"]][[id_row]],
        c_index_cv_sd = fit[["cvsd"]][[id_row]],,
        fit = fit
    )
    return(res)
}

# Fit a standard Cox proportional hazards model
#
# Performs a standard (unpenalized) Cox regression on features and optional
# covariates. This is used when no penalization is desired.
.fit_standard_cox <- function(mat, time, status, col.var, ...){
    # Create survival response object as it is required input format for coxph
    y <- Surv(time, status)
    # Combine survival times, status, features, and optional covariates into one
    # data frame
    df <- data.frame(time = time, status = status, mat, col.var)
    # Fit standard Cox proportional hazards model
    fit <- coxph(y ~ ., data = df)

    # Extract raw coefficients
    coefs <- coef(fit)
    # Filter and normalize coefficients
    coefs <- .filter_and_normalize_coefs(coefs, ...)
    #  Compute risk scores for all samples
    risk_scores <- predict(fit, type = "lp")
    #  Compute apparent concordance index
    c_index <- .compute_c_index(risk_scores, y, fit, penalized = FALSE)

    # Return results as a list
    res <- list(
        coefficients = coefs,
        risk_scores = as.numeric(risk_scores),
        c_index = c_index,
        fit = fit
    )
    return(res)
}

# Filter and normalize coefficients
#
# Removes coefficients with magnitude below a threshold and rescales the
# remaining coefficients.
.filter_and_normalize_coefs <- function(coefs, coef.threshold = 0, ...) {
    # Identify coefficients exceeding the threshold in absolute value
    # to remove very small coefficients so that we can focus on meaningful
    # predictors
    res <- coefs[ abs(coefs) > coef.threshold ]
    # Rescale the remaining coefficients
    if( length(res) > 0L ){
        res <- 2*res / sum(abs(res))
    }
    return(res)
}

# Compute concordance index (C-index) for survival predictions
#
# Evaluates the discriminative ability of a survival model:
# how well the predicted risk scores rank patients according to observed
# survival.
#' @importFrom glmnet Cindex
.compute_c_index <- function(risk_scores, y, fit = NULL, penalized = TRUE){
    res <- NA
    if( !(length(risk_scores) == 0 || all(risk_scores == 0)) ){
        if( penalized ){
            #  Penalized Cox: use glmnet's Cindex function
            res <- Cindex(pred = risk_scores, y)
        } else {
            # Standard Cox: extract C-index from the model summary
            res <- summary(fit)[["concordance"]][[1L]]
        }
    }
    return(res)
}

# Select lambda value from a cross-validated glmnet fit
#
# Allows choosing lambda based on CV results or restricting the model to a
# maximum number of variables.
.lambda_selector <- function(fit, lambda, nvar = NULL, ...) {
    # Handle nvar cutoff
    if( !is.null(nvar) ){
        # Identify lambda values with number of non-zero coefficients <= nvar.
        valid <- which(fit[["glmnet.fit"]][["df"]] <= nvar)
        # If such lambdas exist, choose the lambda for model that has
        # highest number of taxa.
        if( length(valid) > 0L ){
            lambda <- fit[["glmnet.fit"]][max(valid), "lambda"]
        }
    }
    # Handle lambda selection by name
    if( is.character(lambda) ){
        # Convert character string to actual numeric lambda value from fit
        # Why: "lambda.min" or "lambda.1se" are stored inside fit; this
        # retrieves the correct numeric value
        lambda <- fit[[lambda]]
    }
    return(lambda)
}
