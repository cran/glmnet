\name{glmnet.control}
\alias{glmnet.control}
\title{internal glmnet parameters}
\description{
  View and/or change the factory default parameters in glmnet
  }
\usage{
glmnet.control(fdev=1.0e-5, devmax=0.999, eps=1.0e-6, big=9.9e35, mnlam=5, pmin=1.0e-5,
     exmx=250.0,prec=1e-10,mxit=100,factory=FALSE)
}
\arguments{
\item{fdev}{minimum fractional change in deviance for stopping path;
  factory default = 1.0e-5}
\item{devmax}{maximum fraction of explained deviance for stopping path;
  factory default = 0.999}
\item{eps}{minimum value of lambda.min.ratio (see glmnet); factory default= 1.0e-6}
\item{big}{large floating point number; factory default = 9.9e35. Inf in
  definition of upper.limit is set to big}
\item{mnlam}{minimum number of path points (lambda values) allowed;
  factory default = 5}
\item{pmin}{minimum null probability for any class. factory default = 1.0e-5}
\item{exmx}{maximum allowed exponent. factory default = 250.0}
\item{prec}{convergence threshold for multi response bounds adjustment solution. factory default = 1.0e-10}
\item{mxit}{maximum iterations for multiresponse bounds adjustment
  solution.  factory default = 100}
\item{factory}{If \code{TRUE}, reset all the parameters to the factory default;
  default is \code{FALSE}}
}
\details{
  If called with no arguments, \code{glmnet.control()} returns a list
  with the current settings of these parameters. Any arguments included in the
  call sets those parameters to the new values, and then silently
  returns. The values set are persistent for the duration of the R session.
}
\value{
  A list with named elements as in the argument list
}

\author{Jerome Friedman, Trevor Hastie\cr
Maintainer: Trevor Hastie \email{hastie@stanford.edu}}
\seealso{\code{glmnet}}
\examples{
glmnet.control(fdev=0)#continue along path even though not much changes
glmnet.control() # view current settings
glmnet.control(factory=TRUE) # reset all the parameters to their default
}
\keyword{models}
\keyword{regression}
