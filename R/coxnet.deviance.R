#' Compute deviance for Cox model
#'
#' Compute the deviance (-2 log partial likelihood) for Cox model.
#'
#' Computes the deviance for a single set of predictions, or for a matrix
#' of predictions. The user can either supply the predictions
#' directly through the \code{pred} option, or by supplying the \code{x} matrix
#' and \code{beta} coefficients.
#'
#' The function first checks if \code{pred} is passed: if so, it is used as
#' the predictions. If \code{pred} is not passed but \code{x} and \code{beta}
#' are passed, then these values are used to compute the predictions. If
#' neither \code{x} nor \code{beta} are passed, then the predictions are all
#' taken to be 0.
#'
#' Uses the C++ coxdev library for computation, supporting both Breslow
#' and Efron methods for ties, as well as stratified and (start, stop] data.
#'
#' @aliases coxnet.deviance
#' @param pred Fit vector or matrix (usually from glmnet at a particular
#' lambda or a sequence of lambdas).
#' @param y Survival response variable, must be a \code{Surv} or
#' \code{stratifySurv} object.
#' @param x Optional \code{x} matrix, to be supplied if \code{pred = NULL}.
#' @param offset Optional offset vector.
#' @param weights Observation weights (default is all equal to 1).
#' @param std.weights If TRUE (default), observation weights are standardized
#' to sum to 1.
#' @param beta Optional coefficient vector/matrix, to be supplied if
#' \code{pred = NULL}.
#' @param cox.ties Character; the method for handling ties. One of
#'   \code{"breslow"} (the current default) or \code{"efron"}. The default
#'   will change to \code{"efron"} in glmnet 5.1 to match
#'   \code{survival::coxph}.
#'
#' @return A vector of deviances, one for each column of predictions.
#'
#' @examples
#' set.seed(1)
#' eta <- rnorm(10)
#' time <- runif(10, min = 1, max = 10)
#' d <- ifelse(rnorm(10) > 0, 1, 0)
#' y <- survival::Surv(time, d)
#' coxnet.deviance(pred = eta, y = y)
#'
#' # if pred not provided, it is set to zero vector
#' coxnet.deviance(y = y)
#'
#' # example with x and beta
#' x <- matrix(rnorm(10 * 3), nrow = 10)
#' beta <- matrix(1:3, ncol = 1)
#' coxnet.deviance(y = y, x = x, beta = beta)
#'
#' # example with (start, stop] data
#' y2 <- survival::Surv(time, time + runif(10), d)
#' coxnet.deviance(pred = eta, y = y2)
#'
#' # example with strata
#' y2 <- stratifySurv(y, rep(1:2, length.out = 10))
#' coxnet.deviance(pred = eta, y = y2)
#'
#' @seealso \code{coxgrad}
#' @keywords Cox model
#'
#' @export
coxnet.deviance <- function(pred = NULL, y, x = NULL, offset = NULL,
                            weights = NULL, std.weights = TRUE, beta = NULL,
                            cox.ties = c("breslow", "efron")) {
  cox.ties <- match.arg(cox.ties)
  efron <- (cox.ties == "efron")
  y <- response.coxnet(y)
  nobs <- nrow(y)

  # Extract survival times
  if (ncol(y) == 2) {
    start_time <- rep(0.0, nobs)
    stop_time <- as.double(y[, "time"])
    status <- as.integer(y[, "status"])
  } else if (ncol(y) == 3) {
    start_time <- as.double(y[, "start"])
    stop_time <- as.double(y[, "stop"])
    status <- as.integer(y[, "status"])
  } else {
    stop("Response y should have 2 or 3 columns")
  }

  # Extract strata
  if ("strata" %in% names(attributes(y))) {
    strata <- as.integer(as.factor(attr(y, "strata")))
  } else {
    strata <- integer(0)
  }

  # Handle sparse x
  if (!is.null(x) && inherits(x, "sparseMatrix")) {
    if (is.null(beta))
      stop("if x is passed, beta must also be passed")
    pred <- as.matrix(x %*% beta)
    return(coxnet.deviance(pred = pred, y = y, offset = offset,
                           weights = weights, std.weights = std.weights,
                           cox.ties = cox.ties))
  }

  # Resolve pred from inputs
  if (!is.null(pred)) {
    eta <- as.matrix(pred)
  } else if (is.null(x) && is.null(beta)) {
    eta <- matrix(0, nrow = nobs, ncol = 1)
  } else if (!is.null(x) && !is.null(beta)) {
    eta <- as.matrix(x %*% beta)
  } else {
    stop("user must pass either `pred`, or both `x` and `beta`")
  }

  # Add offset
  if (!is.null(offset)) {
    eta <- eta + as.double(offset)
  }

  # Normalize weights
  if (is.null(weights))
    weights <- rep(1, nobs)
  else {
    if (std.weights) weights <- nobs * weights / sum(weights)
    weights <- as.double(weights)
  }

  # Call C++ coxdev
  result <- compute_cox_quantities_exp(
    start = start_time,
    stop = stop_time,
    status = status,
    strata = strata,
    efron = efron,
    eta = eta,
    weights = weights
  )

  result$deviance
}
