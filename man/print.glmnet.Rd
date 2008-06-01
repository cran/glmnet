\name{print.glmnet}
\alias{print.glmnet}
\title{print a glmnet object}
\description{
Print a summary of the glmnet path at each step along the path.
}
\usage{
\method{print}{glmnet}(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{x}{fitted glmnet object}
  \item{digits}{significant digits in printout}
  \item{\dots}{additional print arguments}
}
\details{
  The call that produced the object \code{x} is printed, followed by a three-column
  matrix
  with columns \code{Df}, \code{\%dev} and \code{Lambda}. The \code{Df}
    column is the number of nonzero coefficients (Df is a reasonable
    name only for lasso fits). \code{\%dev} is the percent deviance
    explained (relative to the null deviance).
    }
\value{
The matrix above is silently returned}
\references{Friedman, J., Hastie, T. and Tibshirani, R. (2008)
  \emph{Regularization Paths for Generalized Linear Models via Coordinate
  Descent}}
\author{Jerome Friedman, Trevor Hastie and Rob Tibshirani\cr
  Maintainer: Trevor Hastie <hastie@stanford.edu>}
\seealso{\code{glmnet}, \code{predict} and \code{coef} methods.}
\examples{
x=matrix(rnorm(100*20),100,20)
y=rnorm(100)
fit1=glmnet(x,y)
print(fit1)
}
\keyword{models}
\keyword{regression}

