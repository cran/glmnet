#' print a glmnet object
#'
#' Print a summary of the glmnet path at each step along the path.
#' @details
#' The call that produced the object `x` is printed, followed by a
#' three-column matrix with columns `Df`, `%Dev` and `Lambda`.
#' The `Df` column is the number of nonzero coefficients (Df is a
#' reasonable name only for lasso fits). `%Dev` is the percent deviance
#' explained (relative to the null deviance).  In the case of a 'relaxed' fit,
#' an additional column is inserted, `%Dev R` which gives the percent
#' deviance explained by the relaxed model. For a "bigGlm" model, a simpler
#' summary is printed.
#'
#' @aliases print.glmnet print.relaxed print.bigGlm
#' @param x fitted glmnet object
#' @param digits significant digits in printout
#' @param \dots additional print arguments
#' @return The matrix above is silently returned
#' @seealso \code{glmnet}, \code{predict} and \code{coef} methods.
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008). Regularization Paths for Generalized Linear Models via Coordinate Descent
#' @keywords models regression
#' @examples
#'
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit1 = glmnet(x, y)
#' print(fit1)
#' @method print glmnet
#' @export
print.glmnet <-
    function (x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall: ", deparse(x$call), "\n\n")
     out=data.frame(Df = x$df, `%Dev` = round(x$dev.ratio*100, 2),
                    Lambda = signif(x$lambda, digits),check.names=FALSE,row.names=seq(along=x$df))
class(out)=c("anova",class(out))
     print(out)
}
