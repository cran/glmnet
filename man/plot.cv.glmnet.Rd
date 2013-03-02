\name{plot.cv.glmnet}
\alias{plot.cv.glmnet}
\title{plot the cross-validation curve produced by cv.glmnet}
\description{
Plots the cross-validation curve, and upper and lower standard deviation
curves, as a function of the \code{lambda} values used.}
\usage{
\method{plot}{cv.glmnet}(x, sign.lambda, ...)
}
\arguments{
  \item{x}{fitted \code{"cv.glmnet"} object}
  \item{sign.lambda}{Either plot against \code{log(lambda)} (default) or
    its negative if \code{sign.lambda=-1}.}
  \item{\dots}{Other graphical parameters to plot}
}
\details{
A plot is produced, and nothing is returned.}
\references{Friedman, J., Hastie, T. and Tibshirani, R. (2008)
  \emph{Regularization Paths for Generalized Linear Models via Coordinate
  Descent}}
\author{Jerome Friedman, Trevor Hastie and Rob Tibshirani\cr
  Maintainer: Trevor Hastie <hastie@stanford.edu>}
\seealso{\code{glmnet} and \code{cv.glmnet}.}
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
