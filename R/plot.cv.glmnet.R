#' plot the cross-validation curve produced by cv.glmnet
#'
#' Plots the cross-validation curve, and upper and lower standard deviation
#' curves, as a function of the \code{lambda} values used. If the object has
#' class \code{"cv.relaxed"} a different plot is produced, showing both
#' \code{lambda} and \code{gamma}
#'
#' A plot is produced, and nothing is returned.
#'
#' @aliases plot.cv.glmnet
#' @param x fitted \code{"cv.glmnet"} object
#' @param sign.lambda Either plot against \code{log(lambda)} (default) or its
#' negative if \code{sign.lambda=-1}.
#' @param \dots Other graphical parameters to plot
#' @author Jerome Friedman, Trevor Hastie and Rob Tibshirani\cr Maintainer:
#' Trevor Hastie <hastie@@stanford.edu>
#' @seealso \code{glmnet} and \code{cv.glmnet}.
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008)
#' \emph{Regularization Paths for Generalized Linear Models via Coordinate
#' Descent}
#' @keywords models regression
#' @examples
#'
#' set.seed(1010)
#' n = 1000
#' p = 100
#' nzc = trunc(p/10)
#' x = matrix(rnorm(n * p), n, p)
#' beta = rnorm(nzc)
#' fx = (x[, seq(nzc)] %*% beta)
#' eps = rnorm(n) * 5
#' y = drop(fx + eps)
#' px = exp(fx)
#' px = px/(1 + px)
#' ly = rbinom(n = length(px), prob = px, size = 1)
#' cvob1 = cv.glmnet(x, y)
#' plot(cvob1)
#' title("Gaussian Family", line = 2.5)
#' cvob1r = cv.glmnet(x, y, relax = TRUE)
#' plot(cvob1r)
#' frame()
#' set.seed(1011)
#' par(mfrow = c(2, 2), mar = c(4.5, 4.5, 4, 1))
#' cvob2 = cv.glmnet(x, ly, family = "binomial")
#' plot(cvob2)
#' title("Binomial Family", line = 2.5)
#' ## set.seed(1011)
#' ## cvob3 = cv.glmnet(x, ly, family = "binomial", type = "class")
#' ## plot(cvob3)
#' ## title("Binomial Family", line = 2.5)
#'
#' @method plot cv.glmnet
#' @export
plot.cv.glmnet=function(x,sign.lambda=1,...){
    cvobj=x
    xlab = expression(Log(lambda))
#  xlab="log(Lambda)"
  if(sign.lambda<0)xlab=paste("-",xlab,sep="")
  plot.args=list(x=sign.lambda*log(cvobj$lambda),y=cvobj$cvm,ylim=range(cvobj$cvup,cvobj$cvlo),xlab=xlab,ylab=cvobj$name,type="n")
  new.args=list(...)
  if(length(new.args))plot.args[names(new.args)]=new.args
do.call("plot",plot.args)
    error.bars(sign.lambda*log(cvobj$lambda),cvobj$cvup,cvobj$cvlo,width=0.01,
               col="darkgrey")
#               col="antiquewhite2")

    points(sign.lambda*log(cvobj$lambda),cvobj$cvm,pch=20,
          col="red")
axis(side=3,at=sign.lambda*log(cvobj$lambda),labels=paste(cvobj$nz),tick=FALSE,line=0)
abline(v=sign.lambda*log(cvobj$lambda.min),lty=3)
abline(v=sign.lambda*log(cvobj$lambda.1se),lty=3)
  invisible()
}
