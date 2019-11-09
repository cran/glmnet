#' Extract the deviance from a glmnet object
#'
#' Compute the deviance sequence from the glmnet object
#'
#' A glmnet object has components \code{dev.ratio} and \code{nulldev}.  The
#' former is the fraction of (null) deviance explained. The deviance
#' calculations incorporate weights if present in the model. The deviance is
#' defined to be 2*(loglike_sat - loglike), where loglike_sat is the
#' log-likelihood for the saturated model (a model with a free parameter per
#' observation).  Null deviance is defined to be 2*(loglike_sat
#' -loglike(Null)); The NULL model refers to the intercept model, except for
#' the Cox, where it is the 0 model. Hence dev.ratio=1-deviance/nulldev, and
#' this \code{deviance} method returns (1-dev.ratio)*nulldev.
#'
#' @param object fitted glmnet object
#' @param \dots additional print arguments
#' @return (1-dev.ratio)*nulldev
#' @author Jerome Friedman, Trevor Hastie and Rob Tibshirani\cr Maintainer:
#' Trevor Hastie <hastie@@stanford.edu>
#' @seealso \code{glmnet}, \code{predict}, \code{print}, and \code{coef}
#' methods.
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008)
#' \emph{Regularization Paths for Generalized Linear Models via Coordinate
#' Descent}
#' @keywords models regression
#' @examples
#'
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit1 = glmnet(x, y)
#' deviance(fit1)
#' @method deviance glmnet
#' @export
deviance.glmnet=function(object,...){
dev=object$dev
nulldev=object$nulldev
(1-dev)*nulldev
}
