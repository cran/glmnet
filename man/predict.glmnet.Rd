\name{predict.glmnet}
\alias{coef.glmnet}
\alias{predict.glmnet}
\alias{predict.elnet}
\alias{predict.lognet}
\alias{predict.multnet}
\title{make predictions from a "glmnet" object.}
\description{
Similar to other predict methods, this functions predicts fitted values, logits,
coefficients and more from a fitted \code{"glmnet"} object.
}
\usage{
\method{predict}{glmnet}(object, newx, s = object$lambda,
type=c("link","response","coefficients","class","nonzero"), exact =
FALSE, ...)
\method{coef}{glmnet}(object,s=object$lambda, exact=FALSE, ...)
}
\arguments{
  \item{object}{Fitted \code{"glmnet"} model object.}
  \item{newx}{Matrix of new values for \code{x} at which predictions are
  to be made. Must be a matrix; can be sparse as in \code{Matrix}
  package. This argument is not used for \code{type=c("coefficients","nonzero")}}
  \item{s}{Value(s) of the penalty parameter \code{lambda} at which
    predictions are required. Default is the entire sequence used to
    create the model. }
  \item{type}{Type of prediction required. Type \code{"link"} gives the
    linear predictors for \code{"binomial"} or \code{"multinomial"}
    models; for  \code{"gaussian"} models it gives the fitted
    values. Type \code{"response"} gives the fitted probabilities for
    \code{"binomial"} or \code{"multinomial"}; for \code{"gaussian"}
    type \code{"response"} is equivalent to type \code{"link"}. Type
    \code{"coefficients"} computes the coefficients at the requested
    values for \code{s}.  Note that for
    \code{"binomial"} models, results are returned only for the class
    corresponding to the second level of the factor response.
    Type \code{"class"} applies only to
    \code{"binomial"} or \code{"multinomial"} models, and produces the
    class label corresponding to the maximum probability. Type
    \code{"nonzero"} returns a list of the indices of the nonzero
    coefficients for each value of \code{s}.}
  \item{exact}{By default (\code{exact=FALSE}) the predict function uses linear interpolation
  to make predictions for values of \code{s} that do not coincide with
  those used in the fitting algorithm. Currently \code{exact=TRUE} is
  not implemented, but prints an error message telling the user how to
  achieve the exact predictions. This is done my rerunning the algorithm
with the desired values interspersed (in order) with the values used in
the original fit. This is easily achieved via the R command
\code{lamba=sort(c(object$lambda, new.lambda))}}
  \item{\dots}{Not used. Other arguments to predict. } }
\details{The shape of the objects returned are different for
  \code{"multinomial"} objects. This function actually calls
  \code{NextMethod()},
and the appropriate predict method is invoked for each of the three
model types. \code{coef(...)} is equivalent to \code{predict(type="coefficients",...)}}
\value{The object returned depends on type.}
\references{Friedman, J., Hastie, T. and Tibshirani, R. (2008)
  \emph{Regularization Paths for Generalized Linear Models via Coordinate
  Descent}}
\author{Jerome Friedman, Trevor Hastie and Rob Tibshirani\cr
  Maintainer: Trevor Hastie <hastie@stanford.edu>}
\seealso{\code{glmnet}, and \code{print}, and \code{coef} methods.}
\examples{
x=matrix(rnorm(100*20),100,20)
y=rnorm(100)
g2=sample(1:2,100,replace=TRUE)
g4=sample(1:4,100,replace=TRUE)
fit1=glmnet(x,y)
predict(fit1,newx=x[1:5,],s=c(0.01,0.005))
predict(fit1,type="coef")
fit2=glmnet(x,g2,family="binomial")
predict(fit2,type="response",newx=x[2:5,])
predict(fit2,type="nonzero")
fit3=glmnet(x,g4,family="multinomial")
predict(fit3,newx=x[1:3,],type="response",s=0.01)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
