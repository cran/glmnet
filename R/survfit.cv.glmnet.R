#' Compute a survival curve from a cv.glmnet object
#'
#' Computes the predicted survivor function for a Cox proportional hazards
#' model with elastic net penalty from a cross-validated glmnet model.
#'
#' This function makes it easier to use the results of cross-validation
#' to compute a survival curve.
#'
#' @param formula A class \code{cv.glmnet} object. The object should have
#' been fit with \code{family = "cox"}.
#' @param s Value(s) of the penalty parameter lambda at which predictions
#' are required. Default is the value s="lambda.1se" stored on the CV object.
#' Alternatively s="lambda.min" can be used. If s is numeric, it is taken
#' as the value(s) of lambda to be used.
#' @param ... Other arguments to be passed to \code{survfit.coxnet}.
#'
#' @return If \code{s} is a single value, an object of class "survfitcox"
#' and "survfit" containing one or more survival curves. Otherwise, a list
#' of such objects, one element for each value in \code{s}.
#' Methods defined for survfit objects are print, summary and plot.
#'
#' @examples
#' set.seed(2)
#' nobs <- 100; nvars <- 15
#' xvec <- rnorm(nobs * nvars)
#' x <- matrix(xvec, nrow = nobs)
#' beta <- rnorm(nvars / 3)
#' fx <- x[, seq(nvars / 3)] %*% beta / 3
#' ty <- rexp(nobs, exp(fx))
#' tcens <- rbinom(n = nobs, prob = 0.3, size = 1)
#' y <- survival::Surv(ty, tcens)
#' cvfit <- cv.glmnet(x, y, family = "cox")
#' # default: s = "lambda.1se"
#' survival::survfit(cvfit, x = x, y = y)
#'
#' # s = "lambda.min"
#' survival::survfit(cvfit, s = "lambda.min", x = x, y = y)
#' @importFrom survival survfit
#' @method survfit cv.glmnet
#' @export
survfit.cv.glmnet <- function(formula, s = c("lambda.1se", "lambda.min"), ...) {
  this.call <- match.call()
  object <- formula

  # check that a coxnet model was fit
  if (!("coxnet" %in% class(object$glmnet.fit)))
    stop("survfit only available for Cox models")

  # if user didn't specify s, pick it out from its lambda sequence
  if (is.numeric(s)) lambda <- s
  else
    if (is.character(s)) {
      s <- match.arg(s)
      lambda <- object[[s]]
    }
  else stop("Invalid form for s")
  sf <- survfit.coxnet(object$glmnet.fit, s = lambda, ...)
  sf$call <- this.call
  return(sf)
}
