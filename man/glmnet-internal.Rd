\name{glmnet-internal}
\title{Internal glmnet functions}
\alias{fix.lam}
\alias{error.bars}
\alias{getmin}
\alias{cvlognet}
\alias{cvelnet}
\alias{lambda.interp}
\alias{nonzeroCoef}
\alias{softmax}
\alias{jerr}
\alias{plotCoef}
\alias{zeromat}
\description{Internal glmnet functions}
\usage{
fix.lam(lam)
lambda.interp(lambda,s)
nonzeroCoef(beta, bystep=FALSE)
softmax(x,gap)
}
\author{Trevor Hastie}
\details{These are not intended for use by users. \code{lambda.interp}
  does linear interpolation of the lambdas to obtain a prediction at a
  new point s. \code{softmax} does the classification for multinomial models.
  \code{nonzeroCoef} determines in an efficient manner which variables
  are nonzero in each fit. \code{jerr} prints out error messages from
  the fortran. \code{plotCoef} is called by the \code{plot} method for
  \code{glmnet} objects. 
}
\keyword{internal}
