\name{cv.glmnet}
\alias{cv.glmnet}
\title{Cross-validation for glmnet}
\description{Does k-fold cross-validation for glmnet, produces a plot,
  and returns a value for \code{lambda}}
\usage{
cv.glmnet(x, y, ..., nfolds, foldid, type)
}
\arguments{
  \item{x}{\code{x} matrix as in \code{glmnet}.}
  \item{y}{response \code{y} as in \code{glmnet}.}
  \item{\dots}{Other arguments that can be passed to \code{glmnet}.}
  \item{nfolds}{number of folds - default is 10.}
  \item{foldid}{an optional vector of values between 1 and \code{nfold}
    identifying whhat fold each observation is in. If supplied,
    \code{nfold} can be missing.}
  \item{type}{loss to use for cross-validation. Currently two
    options. The default is \code{type="response"}, which uses
    squared-error for gaussian models, and deviance for logistic
    regression. \code{type="class"} applies to logistic regression only,
    and gives misclassification error.}
}
\details{The function runs \code{glmnet} \code{nfolds}+1 times; the
  first to get the \code{lambda} sequence, and then the remainder to
  compute the fit with each of the folds omitted. The error is
  accumulated, and the average error and standard deviation over the
  folds is computed. This function is a preliminary version, since it
  does not allow the full range of data-types for \code{glmnet} yet.
}
\value{an object of class \code{"cv.glmnet"} is returned, which is a
  list with the ingredients of the cross-validation fit.
  \item{lambda}{the values of \code{lambda} used in the fits.}
  \item{cvm}{The mean cross-validated error - a vector of length
    \code{length(lambda)}.}
  \item{cvsd}{estimate of standard error of \code{svm}.}
  \item{cvup}{upper curve = \code{cvm+cvsd}.}
  \item{cvlo}{lower curve = \code{cvm-cvsd}.}
  \item{nzero}{number of non-zero coefficients at each \code{lambda}.}
  \item{name}{a text string indicating type of measure (for plotting
    purposes).}
  \item{lambda.min}{value of \code{lambda} that gives minimum
    \code{cvm}.}
  \item{lambda.1se}{largest value of \code{lambda} such that error is
    within 1 standard error of the minimum.}
}
\references{Friedman, J., Hastie, T. and Tibshirani, R. (2008)
  \emph{Regularization Paths for Generalized Linear Models via Coordinate
  Descent}\url{http://www-stat.stanford.edu/~hastie/Papers/glmnet.pdf}}
\author{Jerome Friedman, Trevor Hastie and Rob Tibshirani\cr 
Maintainer: Trevor Hastie \email{hastie@stanford.edu}}
\seealso{\code{glmnet} and \code{plot} method for \code{"cv.glmnet"} object.}
\examples{
set.seed(1010)
n=1000;p=100
nzc=trunc(p/10)
x=matrix(rnorm(n*p),n,p)
beta=rnorm(nzc)
fx= (x[,seq(nzc)] \%*\% beta)
eps=rnorm(n)*5
y=drop(fx+eps)
px=exp(fx)
px=px/(1+px)
ly=rbinom(n=length(px),prob=px,size=1)
cvob1=cv.glmnet(x,y)
plot(cvob1)
title("Gaussian Family",line=2.5)
frame()
set.seed(1011)
par(mfrow=c(2,2),mar=c(4.5,4.5,4,1))
cvob2=cv.glmnet(x,ly,family="binomial")
plot(cvob2)
title("Binomial Family",line=2.5)
set.seed(1011)
cvob3=cv.glmnet(x,ly,family="binomial",type="class")
plot(cvob3)
title("Binomial Family",line=2.5)
}
\keyword{models}
\keyword{regression}

