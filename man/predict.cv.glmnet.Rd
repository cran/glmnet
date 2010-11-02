\name{predict.cv.glmnet}
\alias{coef.cv.glmnet}
\alias{predict.cv.glmnet}
\title{make predictions from a "cv.glmnet" object.}
\description{
This function makes predictions from a cross-validated glmnet model,
using the stored \code{"glmnet.fit"} object, and the optimal value
chosen for \code{lambda}.
}
\usage{
\method{predict}{cv.glmnet}(object, newx, s=c("lambda.1se","lambda.min"),...)
\method{coef}{cv.glmnet}(object,s=c("lambda.1se","lambda.min"),...)
}
\arguments{
  \item{object}{Fitted \code{"cv.glmnet"}  object.}
  \item{newx}{Matrix of new values for \code{x} at which predictions are
  to be made. Must be a matrix; can be sparse as in \code{Matrix}
  package. See documentation for \code{predict.glmnet}.}
  \item{s}{Value(s) of the penalty parameter \code{lambda} at which
    predictions are required. Default is the value \code{s="lambda.1se"} stored
    on the CV \code{object}. Alternatively \code{s="lambda.min"} can be
    used. If \code{s} is numeric, it is taken as the value(s) of
    \code{lambda} to be used.}
  \item{\dots}{Not used. Other arguments to predict. } }
\details{This function makes it easier to use the results of
  cross-validation to make a prediction.}
\value{The object returned depends the \dots argument which is passed on
to the \code{predict} method for \code{glmnet} objects.}
\references{Friedman, J., Hastie, T. and Tibshirani, R. (2008)
  \emph{Regularization Paths for Generalized Linear Models via Coordinate
    Descent, Journal of Statistical Software, Vol. 33, Issue 1, Feb 2010}\cr
  \url{http://www-stat.stanford.edu/~hastie/Papers/glmnet.pdf}\cr
  \url{http://www.jstatsoft.org/v33/i01/}
}
\author{Jerome Friedman, Trevor Hastie and Rob Tibshirani\cr
  Maintainer: Trevor Hastie <hastie@stanford.edu>}
\seealso{\code{glmnet}, and \code{print}, and \code{coef} methods, and \code{cv.glmnet}.}
\examples{
x=matrix(rnorm(100*20),100,20)
y=rnorm(100)
cv.fit=cv.glmnet(x,y)
predict(cv.fit,newx=x[1:5,])
coef(cv.fit)
coef(cv.fit,s="lambda.min")
predict(cv.fit,newx=x[1:5,],s=c(0.001,0.002))
}
\keyword{models}
\keyword{regression}
