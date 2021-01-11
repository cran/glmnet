#' print a cross-validated glmnet object
#'
#' Print a summary of the results of cross-validation for a glmnet model.
#'
#' A summary of the cross-validated fit is produced, slightly different for a
#' 'cv.relaxed' object than for a 'cv.glmnet' object.  Note that a 'cv.relaxed'
#' object inherits from class 'cv.glmnet', so by directly invoking
#' \code{print.cv.glmnet(object)} will print the summary as if
#' \code{relax=TRUE} had not been used.
#'
#' @aliases print.cv.glmnet print.cv.relaxed
#' @param x fitted 'cv.glmnet' object
#' @param digits significant digits in printout
#' @param \dots additional print arguments
#' @author Jerome Friedman, Trevor Hastie and Rob Tibshirani\cr Maintainer:
#' Trevor Hastie <hastie@@stanford.edu>
#' @seealso \code{glmnet}, \code{predict} and \code{coef} methods.
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008)
#' \emph{Regularization Paths for Generalized Linear Models via Coordinate
#' Descent}\cr \url{https://arxiv.org/abs/1707.08692}\cr Hastie, T.,
#' Tibshirani, Robert, Tibshirani, Ryan (2019) \emph{Extended Comparisons of
#' Best Subset Selection, Forward Stepwise Selection, and the Lasso}
#' @keywords models regression
#' @examples
#'
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit1 = cv.glmnet(x, y)
#' print(fit1)
#' fit1r = cv.glmnet(x, y, relax = TRUE)
#' print(fit1r)
#' ## print.cv.glmnet(fit1r)  ## CHECK WITH TREVOR
#' @method print cv.glmnet
#' @export
#' @export print.cv.glmnet
print.cv.glmnet <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall: ", deparse(x$call), "\n\n")

    optlams=c(x$lambda.min,x$lambda.1se)
    which=match(optlams,x$lambda)
    mat = with(x, cbind(optlams, which, cvm[which], cvsd[which], nzero[which]))
    dimnames(mat) = list(c("min", "1se"), c("Lambda", "Index","Measure",
                                            "SE", "Nonzero"))
    cat("Measure:", x$name,"\n\n")

    mat=data.frame(mat,check.names=FALSE)
    class(mat)=c("anova",class(mat))
    print(mat,digits=digits)
}

