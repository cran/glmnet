\name{plot.glmnet}
\alias{plot.glmnet}
\alias{plot.multnet}
\alias{plot.mrelnet}
\title{plot coefficients from a "glmnet" object}
\description{
Produces a coefficient profile plot of the coefficient paths for a
fitted \code{"glmnet"} object.
}
\usage{
\method{plot}{glmnet}(x, xvar = c("norm", "lambda", "dev"), label = FALSE, ...)
\method{plot}{multnet}(x, xvar = c("norm", "lambda", "dev"), label = FALSE,type.coef=c("coef","2norm"), ...)
\method{plot}{mrelnet}(x, xvar = c("norm", "lambda", "dev"), label = FALSE,type.coef=c("coef","2norm"), ...)
}
\arguments{
  \item{x}{fitted \code{"glmnet"} model}
  \item{xvar}{What is on the X-axis. \code{"norm"} plots against the
    L1-norm of the coefficients, \code{"lambda"} against the log-lambda
    sequence, and \code{"dev"} against the percent deviance explained.}
  \item{label}{If \code{TRUE}, label the curves with variable sequence
    numbers.}
  \item{type.coef}{If \code{type.coef="2norm"} then a single curve per variable,
    else if \code{type.coef="coef"}, a coefficient plot per response}
  \item{\dots}{Other graphical parameters to plot}
}
\details{
A coefficient profile plot is produced. If \code{x} is a multinomial
model, a coefficient plot is produced for each class.
}
\references{Friedman, J., Hastie, T. and Tibshirani, R. (2008)
  \emph{Regularization Paths for Generalized Linear Models via Coordinate
  Descent}}
\author{Jerome Friedman, Trevor Hastie and Rob Tibshirani\cr
  Maintainer: Trevor Hastie <hastie@stanford.edu>}
\seealso{\code{glmnet}, and \code{print}, \code{predict} and \code{coef} methods.}
\examples{
x=matrix(rnorm(100*20),100,20)
y=rnorm(100)
g2=sample(1:2,100,replace=TRUE)
g4=sample(1:4,100,replace=TRUE)
fit1=glmnet(x,y)
plot(fit1)
plot(fit1,xvar="lambda",label=TRUE)
fit3=glmnet(x,g4,family="multinomial")
plot(fit3,pch=19)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
