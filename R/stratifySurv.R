#' @export
`[.stratifySurv` <- function(x, i, j, drop = FALSE) {
  strata <- attr(x, "strata")
  stop_time <- NULL; start_time <- NULL; ss_match <- NULL
  if ("stop_time" %in% names(attributes(x)))  stop_time  <- attr(x, "stop_time")
  if ("start_time" %in% names(attributes(x))) start_time <- attr(x, "start_time")
  if ("ss_match" %in% names(attributes(x)))   ss_match   <- attr(x, "ss_match")
  obj <- NextMethod(`[`)
  # !missing(i) && missing(j)?
  if (!missing(i) && is.matrix(obj)) {
    attr(obj, "strata") <- strata[i]
    if (!is.null(stop_time))  attr(obj, "stop_time")  <- stop_time[i]
    if (!is.null(start_time)) attr(obj, "start_time") <- start_time[i]
    if (!is.null(ss_match))   attr(obj, "ss_match")   <- ss_match[i]
    class(obj) <- class(x)
  }
  return(obj)
}

#' Add strata to a Surv object
#'
#' Helper function to add strata as an attribute to a Surv object. The
#' output of this function can be used as the response in \code{glmnet()}
#' for fitting stratified Cox models.
#'
#' When fitting a stratified Cox model with \code{glmnet()}, strata should
#' be added to a \code{Surv} response with this helper function. Note that
#' it is not sufficient to add strata as an attribute to the \code{Surv}
#' response manually: if the result does not have class \code{stratifySurv},
#' subsetting of the response will not work properly.
#'
#' @param y A Surv object.
#' @param strata A vector of length equal to the number of observations in
#' y, indicating strata membership. Default is all belong to same strata.
#'
#' @return An object of class \code{stratifySurv} (in addition to all the
#' classes \code{y} belonged to).
#'
#' @examples
#' y <- survival::Surv(1:10, rep(0:1, length.out = 10))
#' strata <- rep(1:3, length.out = 10)
#' y2 <- stratifySurv(y, strata)  # returns stratifySurv object
#'
#' @importFrom survival is.Surv
#' @export
stratifySurv <- function(y, strata = rep(1, length(y))) {
  y <- response.coxnet(y)
  if (length(y) != length(strata))
    stop("y and strata must have the same length (=nobs)")
  attr(y, "strata") <- strata
  y_class <- class(y)
  if (!("stratifySurv" %in% y_class)) class(y) <- c("stratifySurv", y_class)
  return(y)
}
