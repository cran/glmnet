#' Compute gradient for Cox model
#'
#' Compute the gradient of the log partial likelihood at a particular fit for Cox
#' model.
#'
#' Compute a gradient vector at the fitted vector for the log partial likelihood.
#' This is like a residual vector, and useful for manual screening of
#' predictors for \code{glmnet} in applications where \code{p} is very large
#' (as in GWAS).
#'
#' Uses the C++ coxdev library for computation, supporting both Breslow
#' and Efron methods for ties, as well as stratified and (start, stop] data.
#'
#' @aliases coxgrad
#' @param eta Fit vector (usually from glmnet at a particular lambda).
#' @param y Survival response variable, must be a \code{Surv} or
#' \code{stratifySurv} object.
#' @param w Observation weights (default is all equal to 1).
#' @param std.weights If TRUE (default), observation weights are standardized
#' to sum to 1.
#' @param diag.hessian If \code{TRUE}, compute the diagonal of the Hessian
#' of the log partial likelihood as well. Default is \code{FALSE}.
#' @param cox.ties Character; the method for handling ties. One of
#'   \code{"breslow"} (the current default) or \code{"efron"}. The default
#'   will change to \code{"efron"} in glmnet 5.1 to match
#'   \code{survival::coxph}.
#'
#' @return A single gradient vector the same length as \code{eta}. If
#' \code{diag.hessian=TRUE}, the diagonal of the Hessian is
#' included as an attribute "diag_hessian".
#'
#' @examples
#' set.seed(1)
#' eta <- rnorm(10)
#' time <- runif(10, min = 1, max = 10)
#' d <- ifelse(rnorm(10) > 0, 1, 0)
#' y <- survival::Surv(time, d)
#' coxgrad(eta, y)
#'
#' # return diagonal of Hessian as well
#' coxgrad(eta, y, diag.hessian = TRUE)
#'
#' # example with (start, stop] data
#' y2 <- survival::Surv(time, time + runif(10), d)
#' coxgrad(eta, y2)
#'
#' # example with strata
#' y2 <- stratifySurv(y, rep(1:2, length.out = 10))
#' coxgrad(eta, y2)
#'
#' @seealso \code{coxnet.deviance}
#' @keywords Cox model
#'
#' @export
coxgrad <- function(eta, y, w, std.weights = TRUE, diag.hessian = FALSE,
                    cox.ties = c("breslow", "efron")) {
    cox.ties <- match.arg(cox.ties)
    efron <- (cox.ties == "efron")
    nobs <- nrow(y)
    if (missing(w)) w <- rep(1, nobs)
    if (std.weights) w <- w / sum(w)

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

    # Call C++ coxdev
    result <- compute_cox_quantities_exp(
        start = start_time,
        stop = stop_time,
        status = status,
        strata = strata,
        efron = efron,
        eta = as.matrix(as.double(eta)),
        weights = as.double(w)
    )

    grad <- result$gradient
    if (diag.hessian) {
        attr(grad, "diag_hessian") <- result$diag_hessian
    }
    grad
}

#' Helper function for Cox deviance and gradient
#'
#' Helps to find ties in death times of data.
#'
#' @param x Sorted vector of death times.
#' @param index Vector of indices for the death times.
#'
#' @return A list with two arguments.
#' \item{index_first}{A vector of indices for the first observation at each
#' death time as they appear in the sorted list.}
#' \item{index_ties}{If there are no ties at all, this is NULL. If not, this is
#' a list with length equal to the number of unique times with ties. For each
#' time with ties, index_ties gives the indices of the observations with a
#' death at that time.}
#'
#' @examples
#' # Example with no ties
#' glmnet:::fid(c(1, 4, 5, 6), 1:5)
#'
#' # Example with ties
#' glmnet:::fid(c(1, 1, 1, 2, 3, 3, 4, 4, 4), 1:9)
fid <- function(x,index) {
    idup=duplicated(x)
    if(!any(idup)) list(index_first=index,index_ties=NULL)
    else {
        ndup=!idup
        xu=x[ndup]# first death times
        index_first=index[ndup]
        ities=match(x,xu)
        index_ties=split(index,ities)
        nties=sapply(index_ties,length)
        list(index_first=index_first,index_ties=index_ties[nties>1])
    }
}
