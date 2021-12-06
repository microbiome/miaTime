#' \linkS4class{SummarizedExperiment} class object stores time information about
#' each sample in `colData` field. Time information can be categorized by
#' per subject(individual) and then used with the selected time operation.
#' The resulting vector is added as a new column of the`colData` field.
#'
#' @param se : \linkS4class{SummarizedExperiment} class object
#' @param time_field : column vector containing the time information
#' @param subject_field : column vector containing the subject(individual)
#' information
#' @param operation : a function applied to `time_field`
#' @param new_field : character class object refers to the name of a
#' newly added field
#' @param ... for new parameters
#'
#' @examples
#' library(miaTime)
#' data(hitchip1006)
#' se <- hitchip1006
#'
#' hitchipTime <- timePerSubject(se, time_field = "time",
#'     subject_field = "subject", operation ="rev", new_field = "reversed", decreasing = TRUE)
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom S4Vectors DataFrame
#'
#' @return a list object containing each subject individually with the
#' related data provided in the `colData` field
#'
#' @export
timePerSubject <- function(se, time_field, subject_field, operation, new_field, ...){

    colData(se)$time_field <- colData(se)[,time_field]

    colData(se)$subject_field <- colData(se)[, subject_field]

    colData(se) <- colData(se) %>%
        as.data.frame() %>%
        group_by(col = subject_field) %>%
        mutate(new_field = match.fun(operation)(time_field, ...)) %>%
        DataFrame()


    colnames(colData(se))[which(names(colData(se)) == "col")] <- paste(subject_field, "_2", sep = "")

    colnames(colData(se))[which(names(colData(se)) == "new_field")] <- new_field

    drop <- c("time_field","subject_field")
    df <- colData(se)[,!(names(colData(se)) %in% drop)]

    list <- as.list(split(df, df[, subject_field]))

    return(list)

}
