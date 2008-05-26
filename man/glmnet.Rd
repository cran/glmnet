\name{glmnet}
\alias{glmnet}
\title{fit an elasticnet model path}
\description{
Fit a regularization path for the elasticnet at a grid of values for the
regularization parameter lambda. Can deal with all shapes of data,
including very large sparse data matrices. Fits linear,
logistic and multinomial regression models.
}
\usage{
glmnet(x, y, family=c("gaussian","binomial","multinomial"), weights, alpha = 1,
  nlambda = 100, lambda.min = ifelse(nobs<nvars,0.05,0.0001), lambda,
  standardize = TRUE, thresh = 1e-04,  dfmax = nvars + 1,
  pmax = min(dfmax * 1.2, nvars), exclude, penalty.factor = rep(1, nvars),
  maxit=100, HessianExact = FALSE, type = c("covariance", "naive"))
}

\arguments{
  \item{x}{input matrix, of dimension nobs x nvars; each row is an observation vector. Can be in sparse column format (class \code{"dgCMatrix"} as in package \code{Matrix})}
  \item{y}{response variable. Quantitative for \code{family="gaussian"}. For \code{family="binomial"} should be either a factor with two levels, or a two-column matrix of counts or proportions. For \code{family="multinomial"}, can be a \code{nc>=2} level factor, or a matrix with \code{nc} columns of counts or proportions}
  \item{family}{Response type (see above)}
  \item{weights}{observation weights. Can be total counts if responses are proportion matrices. Default is 1 for each observation}
  \item{alpha}{The elasticnet mixing parameter, with
    \eqn{0<\alpha\le 1}. The penalty is defined
    as \deqn{(1-\alpha)||\beta||_2^2+\alpha||\beta||_1.} \code{alpha=1}
    is the lasso penalty; Currently \code{alpha<0.01} not reliable,
    unless you supply your own \code{lambda} sequence}
  \item{nlambda}{The number of \code{lambda} values - default is 100.}
  \item{lambda.min}{Smallest value for \code{lambda}, as a fraction of
    \code{lambda.max}, the (data derived) entry value (i.e. the smallest
  value for which all coefficients are zero). The default depends on the
  sample size \code{nobs} relative to the number of variables
  \code{nvars}. If \code{nobs > nvars}, the default is \code{0.0001},
  close to zero.  If \code{nobs < nvars}, the default is \code{0.05}.
  A very small value of
  \code{lambda.min} will lead to a saturated fit. This is undefined for
  \code{"binomial"} and \code{"multinomial"} models, and \code{glmnet}
  will exit gracefully when the percentage deviance explained is almost
  1.} 
  \item{lambda}{A user supplied \code{lambda} sequence. Typical usage
    is to have the 
    program compute its own \code{lambda} sequence based on
    \code{nlambda} and \code{lambda.min}. Supplying a value of
    \code{lambda} overides this. Use with care - it is better to supply
    a decreasing sequence of \code{lambda} values than a single (small) value.}
  \item{standardize}{Logical flag for variable standardization, prior to
    fitting the model sequence. The coefficients are always returned on
    the original scale. Default is is \code{standardize=TRUE}}
  \item{thresh}{Convergence threshold for coordinate descent. Each inner
  coordinate-descent loop continues until the relative change in any
  coefficient is less than \code{thresh}. Defaults value is \code{1E-4}.}
  \item{dfmax}{Limit the maximum number of variables in the
    model. Useful for very large \code{nvars}, if a partial path is desired.}
  \item{pmax}{Limit the maximum number of variables ever to be nonzero}
  \item{exclude}{Indices of variables to be excluded from the
    model. Default is none. Equivalent to an infinite penalty factor
    (next item).}
  \item{penalty.factor}{Separate penalty factors can be applied to each
    coefficient. This is a number that multiplies \code{lambda} to allow
    differential shrinkage. Can be 0 for some variables, which implies
    no shrinkage, and that variable is always included in the
    model. Default is 1 for all variables (and implictly infinity for
    variables listed in \code{exclude}).}
    \item{maxit}{Maximum number of outer-loop iterations for
      \code{"binomial"} or \code{"multinomial"} families. Default is 100.}
    \item{HessianExact}{Only applies to \code{"binomial"} or
	  \code{"multinomial"} families.
	  If \code{FALSE} (the default), an upper-bound approximation is made to the
    hessian, which is not recalculated at each outer loop.} 
  \item{type}{Two algorithm types are supported for (only)
    \code{family="gaussian"}. The default \code{type="covariance"} saves all
    inner-products ever computed, and can be much faster than
    \code{type="naive"}. The latter can be more efficient for \code{p >>
      N} situations.}
}
\details{
  The sequence of models implied by \code{lambda} is fit by coordinate
  descent. For \code{family="gaussian"} this is the lasso sequence if
  \code{alpha=1}, else it is the elasticnet sequence. For \code{family="binomial"}
  or \code{family="multinomial"}, this is a lasso or elasticnet regularization path
  for fitting the linear logistic or multinomial logistic regression
  paths. Sometimes the sequence is truncated before \code{nlambda}
  values of \code{lambda} have been used, because of instabilities in
  the logistic or multinomial models near a saturated fit. \code{glmnet(...,family="binomial")}
  fits a traditional logistic regression model for the
  log-odds. \code{glmnet(...,family="multinomial")} fits a symmetric multinomial model, where
  each class is represented by a linear model (on the log-scale). The
  penalties take care of redundancies. A two-class \code{"multinomial"} model
  will produce the same fit as the corresponding \code{"binomial"} model,
  except the pair of coefficient matrices will be equal in maxnitude and
  opposite in sign, and half the \code{"binomial"} values. 
}
\value{
An object with S3 class \code{"glmnet","*" }, where \code{"*"} is
\code{"elnet"}, \code{"lognet"} or
\code{"multnet"}  for the three types of models.
  \item{call}{the call that produced this object}
  \item{a0}{Intercept sequence of length \code{length(lambda)}}
  \item{beta}{For \code{"elnet"} and \code{"lognet"} models, a \code{nvars x
      length(lambda)} matrix of coefficients, stored in sparse column
    format (\code{"dgCMatrix"}). For \code{"multnet"}, a list of \code{nc} such
    matrices, one for each class.}
  \item{lambda}{The actual sequence of \code{lambda} values used}
  \item{dev}{The percentage deviance explained (for \code{"elnet"}, this
      is the R-square).}
  \item{nulldev}{Null deviance}
  \item{df}{The number of nonzero coefficients for each value of
    \code{lambda}. For \code{"multnet"}, this is the number of variables
    with a nonzero coefficient for \emph{any} class.}
  \item{dfmat}{For \code{"multnet"} only. A matrix consisting of the
    number of nonzero coefficients per class}
  \item{dim}{dimension of coefficient matrix (ices)}
  \item{npasses}{total passes over the data summed over all lambda
    values}
  \item{jerr}{error flag, for warnings and errors (largely for internal debugging).}
}
\references{Friedman, J., Hastie, T. and Tibshirani, R. (2008)
  \emph{Regularization Paths for Generalized Linear Models via Coordinate
  Descent}}
\author{Jerome Friedman, Trevor Hastie and Rob Tibshirani\cr 
Maintainer: Trevor Hastie <hastie@stanford.edu>}
\seealso{\code{print}, \code{predict} and \code{coef} methods.}
\examples{
x=matrix(rnorm(100*20),100,20)
y=rnorm(100)
g2=sample(1:2,100,replace=TRUE)
g4=sample(1:4,100,replace=TRUE)
fit1=glmnet(x,y)
fit2=glmnet(x,g2,family="binomial")
fit3=glmnet(x,g4,family="multinomial")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}

 
