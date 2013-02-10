\name{glmnet-internal}
\title{Internal glmnet functions}
\alias{auc}
\alias{auc.mat}
\alias{coxnet.deviance}
\alias{cvcompute}
\alias{getcoef}
\alias{getcoef.multinomial}
\alias{response.coxnet}
\alias{fix.lam}
\alias{error.bars}
\alias{getmin}
\alias{elnet}
\alias{mrelnet}
\alias{lognet}
\alias{fishnet}
\alias{coefnorm}
\alias{coxnet}
\alias{cv.lognet}
\alias{cv.elnet}
\alias{cv.multnet}
\alias{cv.mrelnet}
\alias{cv.coxnet}
\alias{cv.fishnet}
\alias{lambda.interp}
\alias{nonzeroCoef}
\alias{rmult}
\alias{glmnet_softmax}
\alias{jerr}
\alias{jerr.elnet}
\alias{jerr.lognet}
\alias{jerr.fishnet}
\alias{jerr.coxnet}
\alias{jerr.mrelnet}
\alias{plotCoef}
\alias{zeromat}
\alias{na.mean}
\description{Internal glmnet functions}
\usage{
fix.lam(lam)
lambda.interp(lambda,s)
nonzeroCoef(beta, bystep=FALSE)
glmnet_softmax(x)
}
\author{Trevor Hastie}
\details{These are not intended for use by users. \code{lambda.interp}
  does linear interpolation of the lambdas to obtain a prediction at a
  new point s. \code{glmnet_softmax} does the classification for multinomial models.
  \code{nonzeroCoef} determines in an efficient manner which variables
  are nonzero in each fit. \code{jerr} prints out error messages from
  the fortran. \code{plotCoef} is called by the \code{plot} method for
  \code{glmnet} objects. 
}
\keyword{internal}
